//
// Copyright (c) 2015 Mahyar Khayatkhoei
// Copyright (c) 2009-2014 Shawn Singh, Glen Berseth, Mubbasir Kapadia, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
//

#include <algorithm>
#include <vector>
#include <util/Geometry.h>
#include <util/Curve.h>
#include <util/Color.h>
#include <util/DrawLib.h>
#include "Globals.h"

using namespace Util;

Curve::Curve(const CurvePoint& startPoint, int curveType) : type(curveType)
{
	controlPoints.push_back(startPoint);
}

Curve::Curve(const std::vector<CurvePoint>& inputPoints, int curveType) : type(curveType)
{
	controlPoints = inputPoints;
	sortControlPoints();
}

// Add one control point to the vector controlPoints
void Curve::addControlPoint(const CurvePoint& inputPoint)
{
	controlPoints.push_back(inputPoint);
	sortControlPoints();
}

// Add a vector of control points to the vector controlPoints
void Curve::addControlPoints(const std::vector<CurvePoint>& inputPoints)
{
	for (int i = 0; i < inputPoints.size(); i++)
		controlPoints.push_back(inputPoints[i]);
	sortControlPoints();
}

// Draw the curve shape on screen, usign window as step size (bigger window: less accurate shape)
void Curve::drawCurve(Color curveColor, float curveThickness, int window)
{
#ifdef ENABLE_GUI
	// Robustness: make sure there is at least two control point: start and end points
	if (!checkRobust()) return;

	// Move on the curve from t=0 to t=finalPoint, using window as step size, and linearly interpolate the curve points
	float time = 0;
	float end = controlPoints[controlPoints.size() - 1].time;
	Point startPoint = controlPoints[0].position;
	Point endPoint;
	float timeWindow = 0.5;
	while (time <= end) {
		if (calculatePoint(endPoint, time)) { 
			DrawLib::drawLine(startPoint, endPoint, curveColor, curveThickness);
			startPoint = endPoint;
		} else {
			std::cerr<<"Failed to find next Point at time " << time << std::endl;
		}
		time += timeWindow;
	}

	return;
#endif
}

bool compareFunction(CurvePoint cp1, CurvePoint cp2) {
	return (cp1.time < cp2.time);
}

// Sort controlPoints vector in ascending order: min-first
void Curve::sortControlPoints()
{

	for(std::vector<CurvePoint>::iterator it = controlPoints.begin(); it < controlPoints.end(); ++it) {
		std::cout << it->time << std::endl;
	}

	std::sort(controlPoints.begin(), controlPoints.end(), compareFunction);
	for(std::vector<CurvePoint>::iterator it = controlPoints.begin(); it < controlPoints.end(); ++it) {
		std::cout << it->time << std::endl;
	}
}

// Calculate the position on curve corresponding to the given time, outputPoint is the resulting position
bool Curve::calculatePoint(Point& outputPoint, float time)
{
	// Robustness: make sure there is at least two control point: start and end points
	if (!checkRobust())
		return false;

	// Define temporary parameters for calculation
	unsigned int nextPoint;
	//float normalTime, intervalTime;

	// Find the current interval in time, supposing that controlPoints is sorted (sorting is done whenever control points are added)
	if (!findTimeInterval(nextPoint, time))
		return false;

	// Calculate position at t = time on curve
	if (type == hermiteCurve)
	{
		outputPoint = useHermiteCurve(nextPoint, time);
	}
	else if (type == catmullCurve)
	{
		outputPoint = useCatmullCurve(nextPoint, time);
	}

	// Return
	return true;
}

// Check Roboustness
bool Curve::checkRobust()
{
	if (controlPoints.size() < 2)
		return false;
	else
		return true;
}

// Find the current time interval (i.e. index of the next control point to follow according to current time)
bool Curve::findTimeInterval(unsigned int& nextPoint, float time)
{

	if (controlPoints[0].time == time) {
		nextPoint = 0;
		return true;
	}
	for (int i = 0; i < controlPoints.size() - 1 ; ++i) {
		if (controlPoints[i].time < time && controlPoints[i+1].time >= time) {
			nextPoint = i + 1;
			return true;
		}
	}

	return false;
}

float h1 (float t) {
	return (2 * t * t * t - 3 * t * t + 1);
}

float h2 (float t) {
	return (-2 * t * t * t + 3 * t * t);
}

float h3 (float t) {
	return (t * t * t - 2 * t * t + t);
}

float h4 (float t) {
	return (t * t * t - t * t);
} 

// Implement Hermite curve
Point Curve::useHermiteCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;
	float normalTime, intervalTime;

	//=========================================================================

	if (nextPoint == 0) {
		newPosition = controlPoints[nextPoint].position;
		return newPosition;
	}
	// Calculate time interval, and normal time required for later curve calculations
	intervalTime = controlPoints[nextPoint].time - controlPoints[nextPoint - 1].time;
	normalTime = (time - controlPoints[nextPoint - 1].time) / intervalTime;
	// Calculate position at t = time on Hermite curve
	Point p0 = controlPoints[nextPoint - 1].position;
	Point p1 = controlPoints[nextPoint].position;
	Vector v0 = controlPoints[nextPoint - 1].tangent;
	Vector v1 = controlPoints[nextPoint].tangent;
	float _h1 = h1(normalTime);
	float _h2 = h2(normalTime);
	float _h3 = h3(normalTime);
	float _h4 = h4(normalTime);
	newPosition.x = p0.x * _h1 + p1.x * _h2 + v0.x * _h3 + v1.x * _h4;
	newPosition.y = p0.y * _h1 + p1.y * _h2 + v0.y * _h3 + v1.y * _h4;
	newPosition.z = p0.z * _h1 + p1.z * _h2 + v0.z * _h3 + v1.z * _h4;
	// Return result
	return newPosition;
}


float c0 (float t, float p0, float p1, float p2, float p3) {
	return p1;
}

float c1 (float t, float p0, float p1, float p2, float p3) {
	return -1 * t * p0 + t * p2;
}

float c2 (float t, float p0, float p1, float p2, float p3) {
	return 2 * t * p0 + (t - 3) * p1 + (3 - 2 * t) * p2 - t * p3;
}
float c3 (float t, float p0, float p1, float p2, float p3) {
	return -1 * t * p0 + (2 - t) * p1 + (t - 2) * p2 + t * p3;
}


float calculateCoordinate(float u, float t, float p0, float p1, float p2, float p3) {
	return c0(t, p0, p1, p2, p3) + c1(t, p0, p1, p2, p3) * u + c2(t, p0, p1, p2, p3) * u * u + c3(t, p0, p1, p2, p3);
}

// Implement Catmull-Rom curve
Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;
	float tension = 0.5;
	if (nextPoint <= 1 || nextPoint >= controlPoints.size() - 1) {
		newPosition = useHermiteCurve(nextPoint, time);
		return newPosition;
	}


	// Calculate time interval, and normal time required for later curve calculations

	//float intervalTime = controlPoints[nextPoint].time - controlPoints[nextPoint - 1].time;
	// Calculate position at t = time on Catmull-Rom curve
	Point p0 = controlPoints[nextPoint - 2].position;
	Point p1 = controlPoints[nextPoint - 1].position;
	Point p2 = controlPoints[nextPoint].position;
	Point p3 = controlPoints[nextPoint + 1].position;

	newPosition.x = calculateCoordinate(time, tension, p0.x, p1.x, p2.x, p3.x);
	newPosition.y = calculateCoordinate(time, tension, p0.y, p1.y, p2.y, p3.y);
	newPosition.z = calculateCoordinate(time, tension, p0.z, p1.z, p2.z, p3.z);
	// Return result
	return newPosition;
}
