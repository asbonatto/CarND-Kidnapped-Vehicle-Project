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
        if (tracer == true){
            cout << "Initialization started...";
        }
        // Creating the normal distribution around GPS coordinates
        default_random_engine gen;
        double std_x = std[0];
        double std_y = std[1];
        double std_theta = std[2];
        
        num_particles = 1;
        
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
            weights.push_back(1.0);
        }
        
        is_initialized = true;
        if (tracer == true){
            cout << "completed." << endl;
        }
        
    }
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
   
   if (tracer == true){
        cout << "Prediction started...";
    }
   
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
        particles[p].theta = theta + noise_theta(gen);
    }
    if (tracer == true){
        cout << "completed." << endl;
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// Finds the predicted measurement that is closest to each observed measurement and assigns the observed measurement to this particular landmark.
    // Assumes same coordinate system
	for (int o = 0; o < observations.size(); o++){
        // Using function defined in helper_functions.h
        double min_dist = std::numeric_limits<double>::infinity();
        for (int p = 0; p < predicted.size(); p++){
            double op_dist = dist(observations[o].x, observations[o].y, predicted[p].x, predicted[p].y);
            if (op_dist < min_dist){
                min_dist = op_dist;
                observations[o].id = predicted[p].id;
            }
        }
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
    
    if (tracer == true){
        cout << "Weight update started...";
    }
	
    double sx = pow(std_landmark[0], 2);
    double sy = pow(std_landmark[1], 2);
    
    double num = 0;
    double den = 2*M_PI*std_landmark[0]*std_landmark[1];
    
    for (int p = 0; p < particles.size(); p++){
        // Particle values for coordinate system transformation
        double px = particles[p].x;
        double py = particles[p].y;
        double ptheta = particles[p].theta;
        
        // Limiting the search to landmarks within sensor range
        vector<LandmarkObs> red_landmarks;
        for (int o = 0; o < observations.size(); o++){
            for (int l = 0; l < map_landmarks.landmark_list.size(); l++){
                LandmarkObs newlndmrk{l, map_landmarks.landmark_list[l].x_f, map_landmarks.landmark_list[l].y_f};
                if (dist(observations[o].x, observations[o].y, map_landmarks.landmark_list[l].x_f, map_landmarks.landmark_list[l].y_f) < sensor_range){
                    red_landmarks.push_back(newlndmrk);
                }
            }
        }
        // Converting the observations to global Csys
        vector<LandmarkObs> obs_in_gcsys;
        for (int o = 0; o < observations.size(); o++){
            
            double x = cos(ptheta)*observations[o].x - sin(ptheta)*observations[o].y + px;
            double y = sin(ptheta)*observations[o].x + cos(ptheta)*observations[o].y + py;
            obs_in_gcsys.push_back(LandmarkObs{ observations[o].id, x, y });
        }
        
        // Mapping the observed with predicted landmarks
        dataAssociation(red_landmarks, obs_in_gcsys);
        
        // Finally updating the weight of the particle
        double weight = 1.0;
        for (int o = 0; o < obs_in_gcsys.size(); o++){
            int id = obs_in_gcsys[o].id;
            double err_x = obs_in_gcsys[o].x - red_landmarks[id].x;
            double err_y = obs_in_gcsys[o].y - red_landmarks[id].y;
            
            // Assuming multivariate Gaussian
            num = exp(-0.5*(pow(err_x, 2)/sx + pow(err_y, 2)/sy));
            weight *= num/den;
        }
        
        particles[p].weight = weight;
        weights[p] = weight;
    }
    if (tracer == true){
        cout << "finished." << endl;
    }
}

void ParticleFilter::resample() {
   	if (tracer == true){
        cout << "Resampling started...";
    }
    default_random_engine gen;
    discrete_distribution<> weights_pmf(weights.begin(), weights.end());
    
    vector<Particle> new_particles;
    for (int p = 0; p < num_particles; p++){
        new_particles.push_back(particles[weights_pmf(gen)]);
    }
    particles = new_particles;
    if (tracer == true){
        cout << "finished." << endl;
    }
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
