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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    if (!is_initialized){
        
        // Creating the normal distribution around GPS coordinates
        default_random_engine gen;
        double std_x = std[0];
        double std_y = std[1];
        double std_theta = std[2];
        
        num_particles = 20;
        
        normal_distribution<double> dist_x(x, std_x);
        normal_distribution<double> dist_y(y, std_y);
        normal_distribution<double> dist_theta(theta, std_theta);
        
        for (int i = 0; i < num_particles; i++){
            // Initialize every particle    
            Particle particle;
            particle.id = i;
            particle.x = dist_x(gen);
            particle.y = dist_y(gen);
            particle.theta = dist_theta(gen);
            particle.weight = 1.0;
            // Put into the vector
            particles.push_back(particle);
        }
        
        is_initialized = true;
    }
    return;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
   
    double EPS = 1E-3;
    
    // Creates the noise distribution
    default_random_engine gen;
    normal_distribution<double> noise_x(0.0, std_pos[0]);
    normal_distribution<double> noise_y(0.0, std_pos[1]);
    normal_distribution<double> noise_theta(0.0, std_pos[2]);
    
    for (int p = 0; p < num_particles; p++){
        // First step is to update a single particle
        double x = particles[p].x;
        double y = particles[p].y;
        double theta = particles[p].theta;
        
        if (fabs(yaw_rate) > EPS) {
            double radius = velocity/yaw_rate;
            double delta_theta = yaw_rate*delta_t;
            x += radius*(sin(theta + delta_theta) - sin(theta));
            y += radius*(cos(theta) - cos(theta + delta_theta));
            theta += delta_theta;
        }
        else{
            x += velocity*delta_t*cos(theta);
            y += velocity*delta_t*sin(theta);
        }
        
        particles[p].x = x + noise_x(gen);
        particles[p].y = y + noise_y(gen);
        particles[p].z = theta + noise_theta(gen);
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// Finds the predicted measurement that is closest to each observed measurement and assigns the observed measurement to this particular landmark.
	for (int o = 0; o < observations.size(); o++){
        // Using function defined in helper_functions.h
        double min_dist = std::numeric_limits<double>::infinity();
        for (int p = 0; p < predicted.size(); p++){
            op_dist = dist(observations[o].x, observations[o].y, predicted[p].x, predicted[p].y);
            if (op_dist < min_dist){
                min_dist = op_dist;
                observations[o].id = predicted[p].id;
            }
        }
        
    }
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
}

void ParticleFilter::resample() {
   	
    default_random_engine gen;
    discrete_distribution<> weights_pmf(weights.begin(), weights.end());
    
    vector<Particle> new_particles;
    for (int p = 0; p < num_particles; p++){
        new_particles.push_back(particles[weights_pmf(gen)]);
    }
    particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
