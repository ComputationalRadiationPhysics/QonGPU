#pragma once

template<typename T>
class Params2D : Params{
	
public:
	
	T get_xmax();
	T get_xmin();
	T get_ymax();
	T get_ymin();
	T get_tmax();
	T get_tmin();
	
private: 
	
	T xmax, xmin;
	T ymax, ymin;
	T tmax, tmin;
};

