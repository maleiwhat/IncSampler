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

	bool InTriangle(Edge* edge, Node2d * target_node ,int &mode);

	Edge* splitTriangle(Edge& edge, Node2d& point);
	Edge* splitOnEdge(Edge& edge, Node2d& point);


	void swapTestEdge(Edge * diagonal);
	


	/*Query Part */
	/// Returns a list of "triangles" (one leading half-edge for each triangle)
	const list<Edge*>& getLeadingEdges() const { return leadingEdges_; }

	/// Returns a list of half-edges (one half-edge for each arc)
	list<Edge*>* getEdges(bool skip_boundary_edges = false) const;

public:

	std::vector<std::array<Handle<Node2d>, 4>> flip_records;
	std::vector<Handle<Node2d>> split_records;



};