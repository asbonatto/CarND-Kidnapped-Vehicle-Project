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
    tracer = false;
    if (!is_initialized){
        
        // Creating the normal distribution around GPS coordinates
        default_random_engine gen;
        if (tracer){
            num_particles = 1;
        }
        else{
            num_particles = 20;
        }
                
        
        
        normal_distribution<double> dist_x(x, std[0]);
        normal_distribution<double> dist_y(y, std[1]);
        normal_distribution<double> dist_theta(theta, std[2]);
        
        if (tracer){
            cout << "INIT ----------------------" << endl;
        }
        
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
            if (tracer){
                cout << "p, x, y, theta, w " << particle.x << ", " << particle.y << ", " << particle.theta << ", " << particle.weight << endl;
            }
        }
        
        is_initialized = true;
        
    }
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
   
    double EPS = 1E-5;
    
    // Creates the noise distribution
    default_random_engine gen;
    normal_distribution<double> noise_x(0.0, std_pos[0]);
    normal_distribution<double> noise_y(0.0, std_pos[1]);
    normal_distribution<double> noise_theta(0.0, std_pos[2]);
    
    if (tracer){
        cout << "PREDICTION ----------------------" << endl;
        cout << "dt, Velocity, Yawrate "  << delta_t << ", " << velocity << ", " << yaw_rate << endl;
    }
    
    for (int p = 0; p < num_particles; p++){
        // First step is to update a single particle
        double x = particles[p].x;
        double y = particles[p].y;
        double theta = particles[p].theta;
        
        if (tracer){
            cout << "#: B(x, y, theta) " << p << ", " << particles[p].x <<", " << particles[p].y << ", "<< particles[p].theta << " | ";
        }
        
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
        
        particles[p].x = x + 1*noise_x(gen);
        particles[p].y = y + 1*noise_y(gen);
        particles[p].theta = theta + 1*noise_theta(gen);
        if (tracer){
            cout << "A(x, y, theta) " << particles[p].x <<", " << particles[p].y << ", "<< particles[p].theta << endl;
        }
    }
    
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// Finds the predicted measurement that is closest to each observed measurement and assigns the observed measurement to this particular landmark.
    // Assumes same coordinate system
    if (tracer){
        cout << "ASSOCIATION ----------------------" << endl;
    }
	for (int o = 0; o < observations.size(); o++){
        // Using function defined in helper_functions.h
        double min_dist = std::numeric_limits<double>::infinity();
        for (int p = 0; p < predicted.size(); p++){
            double op_dist = dist(observations[o].x, observations[o].y, predicted[p].x, predicted[p].y);
            if (op_dist < min_dist){
                min_dist = op_dist;
                observations[o].id = p; // Now the id stores the index directly
            }
        }
        int id = observations[o].id;
        if (tracer){
            cout << "O (# , x, y) " << o << ", " << observations[o].x << ", " << observations[o].y << " | ";
            cout << "L (id, x, y) "   << predicted[id].id << ", " << predicted[id].x << ", " << predicted[id].y << endl;
        }
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	
    double sx = pow(std_landmark[0], 2);
    double sy = pow(std_landmark[1], 2);
    
    double den = 2*M_PI*std_landmark[0]*std_landmark[1];
    
    for (int p = 0; p < particles.size(); p++){
        // Particle values for coordinate system transformation
        double px = particles[p].x;
        double py = particles[p].y;
        double ctheta = cos(particles[p].theta);
        double stheta = sin(particles[p].theta);
        
        // Limiting the search to landmarks within sensor range of particle
        vector<LandmarkObs> red_landmarks;
        for (int l = 0; l < map_landmarks.landmark_list.size(); l++){
            LandmarkObs newlndmrk{map_landmarks.landmark_list[l].id_i
                                , map_landmarks.landmark_list[l].x_f
                                , map_landmarks.landmark_list[l].y_f};
                                
            if (dist(particles[p].x, particles[p].y, map_landmarks.landmark_list[l].x_f, map_landmarks.landmark_list[l].y_f) < sensor_range){
                red_landmarks.push_back(newlndmrk);
            }
        }
        // Converting the observations to global Csys
        vector<LandmarkObs> obs_in_gcsys;
        if (tracer){
            cout << "TRANSFORMATION ----------------------" << endl;
        }
        for (int o = 0; o < observations.size(); o++){
            double x = ctheta*observations[o].x - stheta*observations[o].y + px;
            double y = stheta*observations[o].x + ctheta*observations[o].y + py;
            obs_in_gcsys.push_back(LandmarkObs{observations[o].id, x, y});
            // Debugging the coordinate transformation
            if (tracer){
                cout << "L(x, y) " << observations[o].x << ", " << observations[o].y << " | ";
                cout << "G(x, y) " <<                 x << ", " <<                 y << endl;
            }
        }
        
        // Mapping the observed with predicted landmarks
        dataAssociation(red_landmarks, obs_in_gcsys);
        if (tracer){
            cout << "WEIGHT CALCULATION ----------------------" << endl;
            cout << "stdx, stdy " << std_landmark[0] << ", " << std_landmark[1] << endl;
        }
        
        // Finally updating the weight of the particle
       
        double weight = 1.0;
        for (int o = 0; o < obs_in_gcsys.size(); o++){

        	int id = obs_in_gcsys[o].id;
        	
            double err_x = pow(obs_in_gcsys[o].x - red_landmarks[id].x, 2)/sx;
            double err_y = pow(obs_in_gcsys[o].y - red_landmarks[id].y, 2)/sy;
            if (tracer){
                cout << "Obs (#, id, x, y, prob) " << o << ", " << red_landmarks[id].id  << ", " << err_x << ", " << err_y << ", " << exp(-0.5*(err_x + err_y))/den << endl;
            }
            
            // Assuming multivariate Gaussian
            weight *= exp(-0.5*(err_x + err_y))/den;
        }
        
        particles[p].weight = weight;
        weights[p] = weight;
        
        if (tracer){
            cout << "Particle after update..." << endl;
            cout << "#, x, y, t, w" << p << ", " << particles[p].x << ", " << particles[p].y << ", " << particles[p].theta << ", " << particles[p].weight << endl;
        }
    }
}

void ParticleFilter::resample() {
   	
  default_random_engine gen;
  discrete_distribution<int> distribution(weights.begin(), weights.end());

  vector<Particle> resample_particles;

  for(int p = 0; p < num_particles; p++){
    resample_particles.push_back(particles[distribution(gen)]);
  
  }
  particles = resample_particles;
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
