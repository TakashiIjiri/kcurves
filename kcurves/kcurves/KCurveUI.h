#pragma once

//unmanaged
//singleton


#pragma unmanaged 




#include "kcurves.h"




//stl
#include <vector>

#include <vector>
using namespace std;



#define CIRCLE_R 10




class KCurveUI
{
private:
	KCurveUI();


public:
	inline static KCurveUI *getInst()
	{
		static KCurveUI p;
		return &p;
	}

private:
	bool m_bR,m_bL,m_bM;
	int m_dragPtId;


public:
	vector<EVec2d> m_CPs;
	vector<EVec2d> m_kCurveCP;
	vector<EVec2d> m_curves;

	void LBtnDown (const int x, const int y);
	void RBtnDown (const int x, const int y);
	void MBtnDown (const int x, const int y);
	void LBtnUp   (const int x, const int y);
	void RBtnUp   (const int x, const int y);
	void MBtnUp   (const int x, const int y);
	void MouseMove(const int x, const int y);


private:
	int pickCP(const int x, const int y)
	{
		for( int i=0; i < m_CPs.size(); ++i)
		{
			double d = (m_CPs[i][0] - x) * (m_CPs[i][0] - x) + 
				       (m_CPs[i][1] - y) * (m_CPs[i][1] - y);
			if( d < CIRCLE_R * CIRCLE_R) return i;
		}
		return -1;
	}
};







#pragma managed 







/*
class Vec2
{
public:
	double data[2];

	//constructors
	void Set(const Vec2& src)    { data[0] = src.data[0]; data[1] = src.data[1]; }
	Vec2( double x=0, double y=0){ data[0] = x;           data[1] = y;           }
	Vec2(            const Vec2 &src){ Set(src);               }
	Vec2& operator= (const Vec2 &src){ Set(src); return *this; }


	//operators (�R���X�g���N�^���ĂԂ̂ő��������)
	inline double&       operator[](const int   i )       { return data[i]; }
	inline const double& operator[](const int   i ) const { return data[i]; }
	inline Vec2          operator+ (const Vec2& v ) const { return Vec2( data[0] + v[0], data[1] + v[1] ); }
	inline Vec2          operator- (const Vec2& v ) const { return Vec2( data[0] - v[0], data[1] - v[1] ); }

	inline Vec2& operator+=(const Vec2& v){
		data[0] += v[0];
		data[1] += v[1];
		return *this ;
	}
	inline Vec2& operator-=(const Vec2& v){
		data[0] -= v[0];
		data[1] -= v[1];
		return *this ;
	}

	inline double dot(const Vec2& v){
		return data[0] * v[0] + data[1] * v[1];
	}
	inline double length_sq(){ return data[0] * data[0] + data[1] * data[1] ;}
	inline double length()   { return sqrt( length_sq() ) ; }
	inline void   normalize(){
		double l = length();
		if( l == 0) return;
		data[0] /= l;
		data[1] /= l;
	}
	inline Vec2 normalized(){
		double l = length();
		if( l == 0) return Vec2(0,0);
		return Vec2(data[0] / l, data[1] / l);
	}


	//scholar
	inline        Vec2  operator* (double f) const{ return Vec2( f*data[0], f*data[1]); }
	inline        Vec2  operator/ (double f) const{ return Vec2( data[0]/f, data[1]/f); }
	inline friend Vec2  operator* (double f, const Vec2& vec1){ return vec1*f; }
	inline friend Vec2  operator/ (double f, const Vec2& vec1){ return vec1/f; }
	inline        Vec2& operator*=(double f){ data[0] *= f; data[1] *= f; return *this ;}
	inline        Vec2& operator/=(double f){ data[0] /= f; data[1] /= f; return *this; }
};

*/

