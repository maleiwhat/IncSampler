#pragma once

#include <list>
#include <vector>
#include <array>
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include "..\IncSampler.h"
#include "..\Data\Handle.h"

using namespace std;

/*
1. K - to Blue noise
   找出最小距离峰值 cluster it. 
2. Next K
*/





//struct Point : public virtual HandleId
//{
//	// Data
//	union {
//		struct { real x, y; };
//		struct { real e[2]; };
//	};
//
//
//	real & x(){ return x; }
//	real & y(){ return y; }
//
//
//	// Constructors
//	Point() {};
//	Point(real f) { x = y = f; }
//	Point(real x, real y) { this->x = x; this->y = y; }
//	Point(const real *f) { x = f[0]; y = f[1]; }
//
//	// Operators
//	real  operator[] (int i) const { return e[i]; }
//	real& operator[] (int i) { return e[i]; }
//	const Point& operator+ () const { return *this; }
//	Point  operator- () const {
//		Point p(*this); p.x = -p.x; p.y = -p.y; return p;
//	}
//	Point& operator= (const Point &p) {
//		x = p.x; y = p.y; return *this;
//	}
//	Point& operator+= (const Point &p) {
//		*this = *this + p; return *this;
//	}
//	Point& operator-= (const Point &p) {
//		*this = *this - p; return *this;
//	}
//	Point& operator*= (real f) { *this = *this * f; return *this; }
//	Point& operator/= (real f) { *this = *this / f; return *this; }
//
//	// Friend operators
//	friend bool operator== (const Point &p1, const Point &p2) {
//		if (p1.x != p2.x) return false;
//		if (p1.y != p2.y) return false;
//		return true;
//	}
//	friend bool operator!= (const Point &p1, const Point &p2) {
//		return !(p1 == p2);
//	}
//	friend Point operator* (real f, const Point &p) {
//		Point r(p);
//		r.x *= f;
//		r.y *= f;
//		return r;
//	}
//	friend Point operator* (const Point &p, real f) {
//		return f * p;
//	}
//	friend Point operator/ (const Point &p, real f) {
//		real f_ = 1.f / f;
//		return f_ * p;
//	}
//	friend Point operator+ (const Point &p1, const Point &p2) {
//		Point r;
//		r.x = p1.x + p2.x;
//		r.y = p1.y + p2.y;
//		return r;
//	}
//	friend Point operator- (const Point &p1, const Point &p2) {
//		Point r;
//		r.x = p1.x - p2.x;
//		r.y = p1.y - p2.y;
//		return r;
//	}
//	friend std::istream& operator>> (std::istream &is, Point &p) {
//		for (int i = 0; i < 2; ++i) {
//			std::string s;
//			is >> s;
//			std::replace(s.begin(), s.end(), '(', ' ');
//			std::replace(s.begin(), s.end(), ')', ' ');
//			std::stringstream ss;
//			ss << s;
//			ss >> p.e[i];
//		}
//		return is;
//	}
//	friend std::ostream& operator<< (std::ostream &os, const Point &p) {
//		os << p.x << " " << p.y;
//		return os;
//	}
//
//	// Functions
//	inline bool IsInUnitTorus() const {
//		return (0.f <= x && x < 1.f) && (0.f <= y && y < 1.f);
//	}
//	inline void WrapUnitTorus() {
//		x = (x <  0.f) ? x + 1.f : (x >= 1.f) ? x - 1.f : x;
//		y = (y <  0.f) ? y + 1.f : (y >= 1.f) ? y - 1.f : y;
//	}
//	inline real SquaredDistUnitTorus(const Point &p) const {
//		real x = fabsf(this->x - p.x);
//		real y = fabsf(this->y - p.y);
//		x = (x > .5f) ? 1.f - x : x;
//		y = (y > .5f) ? 1.f - y : y;
//		return x*x + y*y;
//	}
//	inline real DistUnitTorus(const Point &p) const {
//		return sqrtf(SquaredDistUnitTorus(p));
//	}
//};


class Node2d : public virtual HandleId
{
private:
	real x_, y_;
public:

	//	// Constructors
	Node2d() {};
	Node2d(real f) { x_ = y_ = f; }
	//	Point(real x, real y) { this->x = x; this->y = y; }
	//	Point(const real *f) { x = f[0]; y = f[1]; }

	Node2d(real x, real y) { x_ = x; y_ = y; }
	real & x(){ return x_; }
	real & y(){ return y_; }

	Node2d& operator= (const Node2d &p) {
			x_ = p.x_; y_ = p.y_; return *this;
		}
	
	real crossProduct2d(real dx1, real dy1, real dx2, real dy2) {
		return dx1*dy2 - dy1*dx2;
	}
	real scalarProduct2d(real dx1, real dy1, real dx2, real dy2) {
		return dx1*dx2 + dy1*dy2;
	}
};

class Edge 
{
private:
	Handle<Node2d> sourceNode_;
	Edge*  twinEdge_;
	Edge*  nextEdgeInFace_;
	bool isLeadingEdge_;

	//
	// Area
	// Radius

public:
	/// 构造函数
	Edge() {
		sourceNode_ = NULL; twinEdge_ = NULL; nextEdgeInFace_ = NULL;
		isLeadingEdge_ = false;
	}

	/// 析构函数
	~Edge() { if (twinEdge_) twinEdge_->setTwinEdge(NULL); }

	void setSourceNode(Node2d * node) { sourceNode_ = node; }
	void setNextEdgeInFace(Edge * edge) { nextEdgeInFace_ = edge; }
	void setTwinEdge(Edge* edge) { twinEdge_ = edge; }

	/// Sets the edge as a leading edge
	void setAsLeadingEdge(bool val = true) { isLeadingEdge_ = val; }
	/// Checks if an edge is a leading edge
	bool isLeadingEdge() const { return isLeadingEdge_; }

	/// Returns the twin edge
	Edge* getTwinEdge() const { return twinEdge_; };

	/// Returns the next edge in face
	Edge* getNextEdgeInFace() const { return nextEdgeInFace_; }

	/// Retuns the source node
	Node2d* getSourceNode() { return sourceNode_.getPtr(); }

	/// Returns the target node
	Node2d* getTargetNode() { return getNextEdgeInFace()->getSourceNode(); }
};

// Half edge version

class IncrementalBlueNoise
{
protected:
	list<Edge *> leadingEdges_; //  Store the leading edges
	void addLeadingEdge(Edge * edge);
	bool removeLeadingEdgeFromList(Edge* leadingEdge);
public:
	bool inited;
	void cleanAll();
	IncrementalBlueNoise(){ inited = false; }
	~IncrementalBlueNoise(){ cleanAll(); }

	/*Delaunay Part*/

	//Add Bounding Box 
	Edge* initTwoEnclosingTriangles(vector<Node2d*>::iterator first,vector<Node2d*>::iterator last);

	Edge* AddBoundingBox(vector<Node2d*>::iterator first, vector<Node2d*>::iterator last);
	// Add Node 
	void removeTriangle(Edge& edge);

	void createDelaunay(vector<Node2d*>::iterator first, vector<Node2d*>::iterator last);

	bool GlobalInsertNode(Node2d * target_node);

	/*
	   mode : -1 not in the triangle
	          0  in the trianlge
			  1  in the leadingedge
			  2  in the leadingedge->getNext
			  3  in the leadingedge->getNext->getNext
	*/
	bool InTriangle(Edge* edge, Node2d * target_node ,int &mode);

	Edge* splitTriangle(Edge& edge, Node2d& point);
	Edge* splitOnEdge(Edge& edge, Node2d& point);


	void swapTestEdge(Edge * diagonal);
	

	/*Delaunay Part End!!!*/
	
	/**/


	/*View Part*/


	/*Query Part */
	/// Returns a list of "triangles" (one leading half-edge for each triangle)
	const list<Edge*>& getLeadingEdges() const { return leadingEdges_; }

	/// Returns a list of half-edges (one half-edge for each arc)
	list<Edge*>* getEdges(bool skip_boundary_edges = false) const;

public:

	std::vector<std::array<Handle<Node2d>, 4>> flip_records;
	std::vector<Handle<Node2d>> split_records;



};