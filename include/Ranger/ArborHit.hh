#ifndef _POINT_H_
#define _POINT_H_

#include "TVector3.h"
#include <iostream>

class ArborHit
{
	float m_hitTime; 
	float m_depth; 
	int m_hitLayer; 
	TVector3 m_hitPos; 
	int m_subD; 
	int m_stave; 
    float m_energy;

	public:
	
	// ArborHit();
	ArborHit( TVector3 hitPos, int hitLayer, float hitTime, float depth, int stave, int subD , float energy=0.);
	
	void setHit( TVector3 hitPos, int hitLayer, float hitTime, float depth, int stave, int subD ,float energy=0.);
    
    float GetEnergy(){
    return m_energy;
    }

	float GetTime()
	{
		return m_hitTime;
	} 
	int GetLayer()
	{
		return m_hitLayer;
	} 
	TVector3 GetPosition()
	{ 
		return m_hitPos;
	}
	int GetSubD()
	{
		return m_subD; 
	}
	int GetStave()
	{
		return m_stave; 
	}
	float GetDepth()
	{
		return m_depth; 
	}
};

#endif
