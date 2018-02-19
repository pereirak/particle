/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"


using namespace std;

int debug = 0;



void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	LOG("Start Init\n");
	
	default_random_engine gen;
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	num_particles = NUM_PARTICLES;

	for(int i =0; i < num_particles; i++){
		Particle current_particle;
		current_particle.id = i;
		current_particle.x = dist_x(gen);
		current_particle.y = dist_y(gen);
		current_particle.theta = dist_theta(gen);
		current_particle.weight = 1.0;
		  
		particles.push_back(current_particle);
		weights.push_back(current_particle.weight);
	}

	is_initialized=true;
	LOG("Init Success\n");

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	LOG("Start Predict\n");
	
	default_random_engine gen;

	for(int i =0; i < num_particles; i++){

		Particle &p = particles[i];

		double xi = particles[i].x, yi = particles[i].y,  thetai = particles[i].theta;
		double xf, yf, thetaf;

		if(fabs(yaw_rate) < 0.0001){
			xf = xi + (velocity)*delta_t*cos(thetai);
			yf = yi + (velocity)*delta_t*sin(thetai);
			thetaf = thetai;
		}else{
			xf = xi + (velocity/yaw_rate)*(sin(thetai+ yaw_rate*delta_t) - sin(thetai));
			yf = yi + (velocity/yaw_rate)*(cos(thetai)-cos(thetai+ yaw_rate*delta_t));
			thetaf = thetai + yaw_rate*delta_t;
		}
		// normal_distribution<double> dist_x(p.x, std_pos[0]);
		// normal_distribution<double> dist_y(p.y, std_pos[1]);
		// normal_distribution<double> dist_theta(p.theta, std_pos[2]);
		normal_distribution<double> dist_x(xf, std_pos[0]);
		normal_distribution<double> dist_y(yf, std_pos[1]);
		normal_distribution<double> dist_theta(thetaf, std_pos[2]);
		
		particles[i].x =  dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);

	}

	LOG("Predict Successful\n");


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	LOG("Start dataAssociation\n");


	for(int i =0; i < observations.size(); i++){

		LandmarkObs &obs = observations[i];

		double ob_x = observations[i].x, ob_y = observations[i].y;
		LandmarkObs &ob = observations[i];

		int nearest_landmark_id=-1;
		double min_dist = numeric_limits<double>::max();

		for(int j=0; j < predicted.size(); j++){

			double pred_x = predicted[j].x, pred_y = predicted[j].y;

			double dist_n = dist(ob_x,ob_y,pred_x,pred_y);
			
			if(dist_n < min_dist ){

				min_dist = dist_n;
				nearest_landmark_id = predicted[j].id;

			} 

		}

		observations[i].id = nearest_landmark_id;

	}

	LOG("dataAssociation Successful\n");

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html



	LOG("Start updateWeights\n");

	double weights_total = 0;
	for(int i=0; i< num_particles; i++){


		Particle &p = ParticleFilter::particles[i];

		double particle_x = particles[i].x, particle_y = particles[i].y, particle_theta = particles[i].theta;

		vector <LandmarkObs> landmark_obs;
		vector<Map::single_landmark_s> llist = map_landmarks.landmark_list;

		for(int j=0; j < llist.size(); j++){

			Map::single_landmark_s lm = llist[j];

			if ((fabs((particle_x - lm.x_f)) <= sensor_range) && (fabs((particle_y - lm.y_f)) <= sensor_range)) {
        		landmark_obs.push_back(LandmarkObs {lm.id_i, lm.x_f, lm.y_f});
        	}	

		}

		vector <LandmarkObs> predicted;
		vector<LandmarkObs> observations_m;
		vector<LandmarkObs> ob_transformed;


		for (int j= 0; j < observations.size(); j++){

			//xm = xp + cos(thetap)*xc - sin(thetap)*yc;
			//ym = yp + sin(thetap)*xc + cos(thetap)*yc;

			LandmarkObs pred_obs;
			pred_obs.id = -j;
			pred_obs.x = particle_x + cos(particle_theta)*observations[j].x - sin(particle_theta)*observations[j].y;
			pred_obs.y = particle_y + sin(particle_theta)*observations[j].x + cos(particle_theta)*observations[j].y;

			ob_transformed.push_back(pred_obs);
		}
			
		dataAssociation(landmark_obs, ob_transformed);


		particles[i].weight=1.0;

		double sig_x = std_landmark[0], sig_y = std_landmark[1];
		double norm = (1.0/(2.0 * M_PI * sig_x * sig_y));

		vector<int> p_association; 
        vector<double> p_sense_x;
        vector<double> p_sense_y;


        // Calculate Weight of particle
		for(int j=0; j < ob_transformed.size(); j++){

			LandmarkObs landmark_nearest;

			for(int k=0; k < landmark_obs.size();k++){

				if(landmark_obs[k].id == ob_transformed[j].id){

					landmark_nearest.id =  landmark_obs[k].id;
					landmark_nearest.x = landmark_obs[k].x;
					landmark_nearest.y = landmark_obs[k].y;
					
					p_association.push_back(landmark_nearest.id); 
                    p_sense_x.push_back(landmark_nearest.x );
                    p_sense_y.push_back(landmark_nearest.y);

					double x_diff = landmark_nearest.x - ob_transformed[j].x;
					double y_diff = landmark_nearest.y - ob_transformed[j].y;
					double exponent= (x_diff*x_diff)/(2 * sig_x*sig_x) + (y_diff*y_diff)/(2 * sig_y*sig_y);
					double weight = norm*exp(-exponent);

					particles[i].weight*=weight;
				}



			}

		}
		//SetAssociations(particles[i],p_association,p_sense_x,p_sense_y);
		weights_total += particles[i].weight;
	}
	

	//Normalize weights
	for(int i=0; i< ParticleFilter::num_particles; i++){
		particles[i].weight /= weights_total;
		ParticleFilter::weights[i] = particles[i].weight;
	}

	LOG("Weight update Successful\n");

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	LOG("Start Resample\n");

	vector<Particle> new_particles;
	random_device seed;
	mt19937 random_generator(seed());

	discrete_distribution<> sample(weights.begin(), weights.end());

	for (int i = 0; i < num_particles; i++)
	{
		new_particles.push_back(particles[sample(random_generator)]);
	}
	particles = new_particles;


	LOG("Resample Successful\n");
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
