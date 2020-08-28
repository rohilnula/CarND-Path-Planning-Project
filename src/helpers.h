#ifndef HELPERS_H
#define HELPERS_H

#include <math.h>
#include <string>
#include <vector>

#include "spline.h"

// for convenience
using std::string;
using std::vector;
using std::pair;

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
//   else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

//
// Helper functions related to waypoints and converting from XY to Frenet
//   or vice versa
//

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Calculate distance between two points
double distance(double x1, double y1, double x2, double y2) {
  return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

// Calculate closest waypoint to current x, y position
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, 
                    const vector<double> &maps_y) {
  double closestLen = 100000; //large number
  int closestWaypoint = 0;

  for (int i = 0; i < maps_x.size(); ++i) {
    double map_x = maps_x[i];
    double map_y = maps_y[i];
    double dist = distance(x,y,map_x,map_y);
    if (dist < closestLen) {
      closestLen = dist;
      closestWaypoint = i;
    }
  }

  return closestWaypoint;
}

// Returns next waypoint of the closest waypoint
int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, 
                 const vector<double> &maps_y) {
  int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

  double map_x = maps_x[closestWaypoint];
  double map_y = maps_y[closestWaypoint];

  double heading = atan2((map_y-y),(map_x-x));

  double angle = fabs(theta-heading);
  angle = std::min(2*pi() - angle, angle);

  if (angle > pi()/2) {
    ++closestWaypoint;
    if (closestWaypoint == maps_x.size()) {
      closestWaypoint = 0;
    }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, 
                         const vector<double> &maps_x, 
                         const vector<double> &maps_y) {
  int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

  int prev_wp;
  prev_wp = next_wp-1;
  if (next_wp == 0) {
    prev_wp  = maps_x.size()-1;
  }

  double n_x = maps_x[next_wp]-maps_x[prev_wp];
  double n_y = maps_y[next_wp]-maps_y[prev_wp];
  double x_x = x - maps_x[prev_wp];
  double x_y = y - maps_y[prev_wp];

  // find the projection of x onto n
  double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
  double proj_x = proj_norm*n_x;
  double proj_y = proj_norm*n_y;

  double frenet_d = distance(x_x,x_y,proj_x,proj_y);

  //see if d value is positive or negative by comparing it to a center point
  double center_x = 1000-maps_x[prev_wp];
  double center_y = 2000-maps_y[prev_wp];
  double centerToPos = distance(center_x,center_y,x_x,x_y);
  double centerToRef = distance(center_x,center_y,proj_x,proj_y);

  if (centerToPos <= centerToRef) {
    frenet_d *= -1;
  }

  // calculate s value
  double frenet_s = 0;
  for (int i = 0; i < prev_wp; ++i) {
    frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
  }

  frenet_s += distance(0,0,proj_x,proj_y);

  return {frenet_s,frenet_d};
}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, 
                     const vector<double> &maps_x, 
                     const vector<double> &maps_y) {
  int prev_wp = -1;

  while (s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1))) {
    ++prev_wp;
  }

  int wp2 = (prev_wp+1)%maps_x.size();

  double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),
                         (maps_x[wp2]-maps_x[prev_wp]));
  // the x,y,s along the segment
  double seg_s = (s-maps_s[prev_wp]);

  double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
  double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

  double perp_heading = heading-pi()/2;

  double x = seg_x + d*cos(perp_heading);
  double y = seg_y + d*sin(perp_heading);

  return {x,y};
}

//Check if any obstacles in front, left and right
pair < vector<bool>, vector<double> > check_if_vehicles_around(vector<vector<double>> sensor_fusion, int prev_size, double car_s, int car_lane) {
  
  bool obs_ahead = false;
  bool obs_left = false;
  bool obs_right = false;
  double distance_to_collision = 1000000;
  double obstacle_velocity = 0.0;
  pair < vector<bool>, vector<double> > result;
  
  for (int i = 0; i < sensor_fusion.size(); i++) {
    
    // calculate lane of obstacle fro sensor data
    int obs_lane = 0;
    if (sensor_fusion[i][6] > 0 && sensor_fusion[i][6] <= 4) {
      obs_lane = 0;
    }
    else if (sensor_fusion[i][6] > 4 && sensor_fusion[i][6] <= 8) {
      obs_lane = 1;
    }
    else {
      obs_lane = 2;
    }
    // (int) (sensor_fusion[i][6] / 4);
    
    // Get position and velocity of the obstacle
    double obs_vx = sensor_fusion[i][3];
    double obs_vy = sensor_fusion[i][4];
    double obs_s = sensor_fusion[i][5];
    
    // Calculate obstacle velocity and its s position by the time our car has the previous trajectory planned
    double obs_speed = sqrt(obs_vx * obs_vx + obs_vy * obs_vy);
    obs_s += prev_size * 0.02 * obs_speed;
    
    // Mark if obstacle vehicle ahead
    if ( (car_lane == obs_lane) && (obs_s > car_s) && (obs_s - car_s < 30) ) {
      obs_ahead = true;
      if (obs_s - car_s < distance_to_collision) {
        distance_to_collision = obs_s - car_s;
        obstacle_velocity = obs_speed;
      }
    }
    
    // Mark if obstcle vehicle to the left
    if ( (car_lane > 0) && (car_lane - 1 == obs_lane) && std::abs(obs_s - car_s) < 20 ) {
      obs_left = true;
    }
    
    // Mark if obstacle vehicle to the right
    if ( (car_lane < 2) && (car_lane + 1 == obs_lane) && std::abs(obs_s - car_s) < 20 ) {
      obs_right = true;
    }
    
  }
  
  if ( car_lane == 0 ) {
    obs_left = true;
  }
  
  if ( car_lane == 2 ) {
    obs_right = true; 
  }
  
  result.first = {obs_ahead, obs_left, obs_right};
  result.second = {distance_to_collision, obstacle_velocity};
  return result;
}

// Try get to max velocity
vector<double> try_get_to_speed(bool obs_ahead, bool obs_left, bool obs_right, int car_lane, double car_vel, double MAX_VEL, double MAX_ACC, double distance_to_collision, double obstacle_velocity) {
  
  double new_lane = car_lane;
  double new_vel = car_vel;
  
  if (!obs_ahead && car_vel >= MAX_VEL) {
    new_lane = car_lane;
    new_vel = car_vel;
    std::cout << "Maintain Lane" << std::endl;
  }
  else if (!obs_ahead && car_vel < MAX_VEL) {
    new_vel = car_vel + MAX_ACC;
    std::cout << "Increase Speed" << std::endl;
  }
  else {
    if (!obs_left) {
      new_lane = car_lane - 1;
      new_vel = car_vel + MAX_ACC * 2;
      std::cout << "<===========================================Move Left" << std::endl;
    }
    else if (!obs_right) {
      new_lane = car_lane + 1;
      new_vel = car_vel + MAX_ACC * 2;
      std::cout << "Move Right============================================>" << std::endl;
    }
    else {
      new_vel = car_vel - MAX_ACC;
      std::cout << "Reduce Speed" << std::endl;
      if (distance_to_collision < 10 and obstacle_velocity < car_vel) {
        new_vel = car_vel - MAX_ACC;
        std::cout << "Reduce Speed2" << std::endl;
      }
    }
  }
  
  if ( new_vel > MAX_VEL ) {
    new_vel = MAX_VEL;
  }
  
  return {new_lane, new_vel};
}

// Calculate Next Lane and velocity
vector<double> get_resultant_speed_and_lane(bool obs_ahead, bool obs_left, bool obs_right, int car_lane, double car_vel, double MAX_VEL, double MAX_ACC, double distance_to_collision, double obstacle_velocity) {
  
  double new_lane = car_lane;
  double new_vel = car_vel;
  
  auto new_lane_vel = try_get_to_speed(obs_ahead, obs_left, obs_right, car_lane, car_vel, MAX_VEL, MAX_ACC, distance_to_collision, obstacle_velocity);
  new_lane = new_lane_vel[0];
  new_vel = new_lane_vel[1];
  
  return {new_lane, new_vel};
}


// Generate Trajectory
vector<vector<double>> generate_smooth_trajectory(int prev_size, double car_lane, double car_vel, 
                                                                        double car_x, double car_y, double car_yaw, 
                                                                        vector<double> previous_path_x, 
                                                                        vector<double> previous_path_y, double car_s, 
                                                                        const vector<double> &map_waypoints_s, 
                                                                        const vector<double> &map_waypoints_x, 
                                                                        const vector<double> &map_waypoints_y) {

  // Vector of x and y co-oordinates placed 30 meters apart
  vector<double> ptsx;
  vector<double> ptsy;
          
  // reference x, y and yaw. They are either the car's start state or the state at the previous paths end point
  double ref_x = car_x;
  double ref_y = car_y;
  double ref_yaw = deg2rad(car_yaw);
          
  // If the size of previous points is less than 2, use the the car's starting position.
  if ( prev_size < 2 ) {
    // Use two points that make the path tangent to the car
    double prev_car_x = car_x - cos(car_yaw);
    double prev_car_y = car_y - sin(car_yaw);
            
    ptsx.push_back(prev_car_x);
    ptsx.push_back(car_x);
            
    ptsy.push_back(prev_car_y);
    ptsy.push_back(car_y);
  } 
  // Use previous path's end point as starting referene
  else {
    // Redefine previous state as previous path endpoint
    ref_x = previous_path_x[prev_size - 1];
    ref_y = previous_path_y[prev_size - 1];
            
    double ref_x_prev = previous_path_x[prev_size - 2];
    double ref_y_prev = previous_path_y[prev_size - 2];
            
    ref_yaw = atan2(ref_y-ref_y_prev, ref_x-ref_x_prev);
            
    ptsx.push_back(ref_x_prev);
    ptsx.push_back(ref_x);
            
    ptsy.push_back(ref_y_prev);
    ptsy.push_back(ref_y);
  }
          
  // In Frenet add points that are 30m spaced from the starting reference
  vector<double> next_wp0 = getXY(car_s + 30, 2 + 4*car_lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
  vector<double> next_wp1 = getXY(car_s + 60, 2 + 4*car_lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
  vector<double> next_wp2 = getXY(car_s + 90, 2 + 4*car_lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
          
  ptsx.push_back(next_wp0[0]);
  ptsx.push_back(next_wp1[0]);
  ptsx.push_back(next_wp2[0]);

  ptsy.push_back(next_wp0[1]);
  ptsy.push_back(next_wp1[1]);
  ptsy.push_back(next_wp2[1]);
          
  for ( int i = 0; i < ptsx.size(); i++ ) {
    // shift car reference angle to 0 degrees
    double shift_x = ptsx[i] - ref_x;
    double shift_y = ptsy[i] - ref_y;

    ptsx[i] = shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw);
    ptsy[i] = shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw);
  }
          
  // Create the spline.
  tk::spline s;
          
  // Set (x, y) point to spline
  s.set_points(ptsx, ptsy);
          
  // Next set of 50 points
  vector<double> next_x_vals;
  vector<double> next_y_vals;
          
  // All path points from last time
  for ( int i = 0; i < prev_size; i++ ) {
    next_x_vals.push_back(previous_path_x[i]);
    next_y_vals.push_back(previous_path_y[i]);
  }
          
  // Calculate distance to a point where x component is 30m ahead
  double target_x = 30.0;
  double target_y = s(target_x);
  double target_dist = sqrt(target_x*target_x + target_y*target_y);

  double x_add_on = 0;

  // Fill up the next points. Need 50 points as next positions
  for( int i = 1; i < 50 - prev_size; i++ ) {

    double N = target_dist/(0.02*car_vel/2.24);
    double x_point = x_add_on + target_x/N;
    double y_point = s(x_point);

    x_add_on = x_point;

    double x_ref = x_point;
    double y_ref = y_point;

    x_point = x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw);
    y_point = x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw);

    x_point += ref_x;
    y_point += ref_y;

    next_x_vals.push_back(x_point);
    next_y_vals.push_back(y_point);
  }
  return { next_x_vals, next_y_vals };
}
#endif  // HELPERS_H