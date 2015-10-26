#include "IncrementalBlueNoise.h"

static Edge * getLeadingEdgeInTriangle(Edge * e)
{
	Edge  * edge = e;

	if (!edge->isLeadingEdge())
	{
		edge = edge->getNextEdgeInFace();
		if (!edge->isLeadingEdge())
			edge = edge->getNextEdgeInFace();
	}

	if (!edge->isLeadingEdge())
		return NULL;

	return edge;
}

static void getLimits(vector<Node2d * >::iterator first,
	vector<Node2d * >::iterator last,
	double & xmin, double & ymin,
	double & xmax, double & ymax)

{
	xmin = ymin = 1.0e+30;
	xmax = ymax = -1.0e+30;

	vector<Node2d*>::iterator it = first;

	for (; it != last; ++it)
	{
		xmin = min(xmin, (*it)->x());
		ymin = min(ymin, (*it)->y());
		xmax = max(xmax, (*it)->x());
		ymax = max(ymax, (*it)->y());
	}
}

static inline real distance_pow2(Handle<Node2d>   n1, Handle<Node2d>   n2)
{
	return (n1->x() - n2->x()) * (n1->x() - n2->x()) + (n1->y() - n2->y()) * (n1->y() - n2->y());
}
static inline real distance_pow2(Node2d *  n1, Node2d * n2)
{
	return (n1->x() - n2->x()) * (n1->x() - n2->x()) + (n1->y() - n2->y()) * (n1->y() - n2->y());
}


// 如果pa pb pc 是按照逆时针排列则返回正值，否则返回负值。 
// 如果在一条直线上，返回零。 
// 这个函数也可以返回这三个点定义的三角形的面积的两倍的近似。
static inline real CounterClockWise(Node2d * na, Node2d * nb, Node2d * nc)
{
	//return ((pb->x - pa->x)*(pc->y - pb->y) - (pc->x - pb->x)*(pb->y - pa->y));
	return ((nb->x() - na->x())*(nc->y() - nb->y()) - (nc->x() - nb->x())*(nb->y() - na->y()));
}


/*
pb
/\
/  \
/     \
/       \
pa ---------   pc
\       /
\    /
\  /
\/
pd
*/



static real InCircle(Node2d * pa, Node2d * pb, Node2d * pc, Node2d * pd)
{
	real det = 0;

	real adx = pa->x() - pd->x();
	real ady = pa->y() - pd->y();

	real bdx = pb->x() - pd->x();
	real bdy = pb->y() - pd->y();

	real cdx = pc->x() - pd->x();
	real cdy = pc->y() - pd->y();

	real bdx_cdy = bdx * cdy;
	real cdx_bdy = cdx * bdy;

	real alift = adx * adx + ady * ady;

	real cdx_ady = cdx * ady;
	real adx_cdy = adx * cdy;

	real blift = bdx * bdx + bdy * bdy;

	real adx_bdy = adx * bdy;
	real bdx_ady = bdx * ady;

	real clift = cdx * cdx + cdy * cdy;

	det = alift * (bdx_cdy - cdx_bdy)
		+ blift * (cdx_ady - adx_cdy)
		+ clift * (adx_bdy - bdx_ady);

	return -det;
}


static real InCircle(Handle<Node2d> pa, Handle<Node2d> pb, Handle<Node2d> pc, Handle<Node2d> pd)
{
	real det = 0;

	real adx = pa->x() - pd->x();
	real ady = pa->y() - pd->y();

	real bdx = pb->x() - pd->x();
	real bdy = pb->y() - pd->y();

	real cdx = pc->x() - pd->x();
	real cdy = pc->y() - pd->y();

	real bdx_cdy = bdx * cdy;
	real cdx_bdy = cdx * bdy;

	real alift = adx * adx + ady * ady;

	real cdx_ady = cdx * ady;
	real adx_cdy = adx * cdy;

	real blift = bdx * bdx + bdy * bdy;

	real adx_bdy = adx * bdy;
	real bdx_ady = bdx * ady;

	real clift = cdx * cdx + cdy * cdy;

	det = alift * (bdx_cdy - cdx_bdy)
		+ blift * (cdx_ady - adx_cdy)
		+ clift * (adx_bdy - bdx_ady);

	return -det;
}


void IncrementalBlueNoise::addLeadingEdge(Edge * edge)
{
	if (edge == NULL) return;
	edge->setAsLeadingEdge();
	leadingEdges_.push_front(edge);
};

//--------------------------------------------------------------------------------------------------
bool IncrementalBlueNoise::removeLeadingEdgeFromList(Edge* leadingEdge) {

	// Remove the edge from the list of leading edges,
	// but don't delete it.
	// Also set flag for leading edge to false.
	// Must search from start of list. Since edges are added to the
	// start of the list during triangulation, this operation will
	// normally be fast (when used in the triangulation algorithm)

	if (leadingEdge == NULL) return false;

#ifdef _DEBUG	
	int i = 0;
#endif 
	list<Edge*>::iterator it;
	for (it = leadingEdges_.begin(); it != leadingEdges_.end(); ++it) {

#ifdef _DEBUG	
		++i;
#endif 
		Edge* edge = *it;
		if (edge == leadingEdge) {

			edge->setAsLeadingEdge(false);
			it = leadingEdges_.erase(it);

			break;
		}
	}

	if (it == leadingEdges_.end())
		return false;

	return true;
}
void IncrementalBlueNoise::cleanAll()
{
	list<Edge*>::const_iterator it;
	for (it = leadingEdges_.begin(); it != leadingEdges_.end(); ++it)
	{
		Edge* e1 = *it;
		Edge* e2 = e1->getNextEdgeInFace();
		Edge* e3 = e2->getNextEdgeInFace();

		delete e1;
		delete e2;
		delete e3;
	}

	leadingEdges_.clear();
	inited = false;
}

//--------------------------------------------------------------------------------------------------
Edge* IncrementalBlueNoise::initTwoEnclosingTriangles(vector<Node2d*>::iterator first,
	vector<Node2d*>::iterator last)
{
	double xmin, ymin, xmax, ymax = 0.0;
	getLimits(first, last, xmin, ymin, xmax, ymax);

	// Add 10% of range:  
	double fac = 10.0;
	double dx = (xmax - xmin) / fac;
	double dy = (ymax - ymin) / fac;

	Node2d* n1 = new Node2d(xmin - dx, ymin - dy);
	Node2d* n2 = new Node2d(xmax + dx, ymin - dy);
	Node2d* n3 = new Node2d(xmax + dx, ymax + dy);
	Node2d* n4 = new Node2d(xmin - dx, ymax + dy);

	// diagonal
	Edge* e1d = new Edge; // lower
	Edge* e2d = new Edge; // upper, the twin edge

	// lower triangle
	Edge* e11 = new Edge;
	Edge* e12 = new Edge;

	// upper triangle
	Edge* e21 = new Edge; // upper upper
	Edge* e22 = new Edge;

	// lower triangle
	e1d->setSourceNode(n3);
	e1d->setNextEdgeInFace(e11);
	e1d->setTwinEdge(e2d);
	e1d->setAsLeadingEdge();
	addLeadingEdge(e1d);

	e11->setSourceNode(n1);
	e11->setNextEdgeInFace(e12);

	e12->setSourceNode(n2);
	e12->setNextEdgeInFace(e1d);

	// upper triangle
	e2d->setSourceNode(n1);
	e2d->setNextEdgeInFace(e21);
	e2d->setTwinEdge(e1d);
	e2d->setAsLeadingEdge();
	addLeadingEdge(e2d);

	e21->setSourceNode(n3);
	e21->setNextEdgeInFace(e22);

	e22->setSourceNode(n4);
	e22->setNextEdgeInFace(e2d);
	
	inited = true;
	
	return e11;
}

Edge* IncrementalBlueNoise::AddBoundingBox(vector<Node2d*>::iterator first,
	vector<Node2d*>::iterator last)
{
	double xmin, ymin, xmax, ymax = 0.0;
	getLimits(first, last, xmin, ymin, xmax, ymax);

	real  max = 0;
	max = xmax > ymax ? xmax : ymax;
	const int size = 1;

	Node2d* n1 = new Node2d(-size * max, 0);
	Node2d* n2 = new Node2d(2*size * max, 0);
	Node2d* n3 = new Node2d(0.5, size * max);
	
	// diagonal
	Edge* e1 = new Edge; // lower
	Edge* e2 = new Edge; // upper, the twin edge
	Edge* e3 = new Edge;

	e1->setSourceNode(n1);
	e2->setSourceNode(n2);
	e3->setSourceNode(n3);

	e1->setNextEdgeInFace(e2);
	e2->setNextEdgeInFace(e3);
	e3->setNextEdgeInFace(e1);

	e1->setAsLeadingEdge();
	addLeadingEdge(e1);

	inited = true;

	return e1;
}




//--------------------------------------------------------------------------------------------------
void IncrementalBlueNoise::removeTriangle(Edge& edge) {

	Edge* e1 = getLeadingEdgeInTriangle(&edge);
#ifdef DEBUG_HE
	if (e1 == NULL)
		errorAndExit("Triangulation::removeTriangle: could not find leading edge");
#endif

	removeLeadingEdgeFromList(e1);
	// cout << "No leading edges = " << leadingEdges_.size() << endl;  
	// Remove the triangle
	Edge* e2 = e1->getNextEdgeInFace();
	Edge* e3 = e2->getNextEdgeInFace();

	if (e1->getTwinEdge())
		e1->getTwinEdge()->setTwinEdge(NULL);
	if (e2->getTwinEdge())
		e2->getTwinEdge()->setTwinEdge(NULL);
	if (e3->getTwinEdge())
		e3->getTwinEdge()->setTwinEdge(NULL);

	delete e1;
	delete e2;
	delete e3;
}

void IncrementalBlueNoise::createDelaunay(vector<Node2d*>::iterator first, vector<Node2d*>::iterator last)
{
	Edge * bedge = initTwoEnclosingTriangles(first, last);
	//Edge * bedge = AddBoundingBox(first, last);

	vector<Node2d*>::iterator it;
	bool status;
	for (it = first; it != last; ++it) {
		status = GlobalInsertNode(*it);
	}


	inited = true;
}


//--------------------------------------------------------------------------------------------------
list<Edge*>* IncrementalBlueNoise::getEdges(bool skip_boundary_edges) const {

	// collect all arcs (one half edge for each arc)
	// (boundary edges are also collected).

	list<Edge*>::const_iterator it;
	list<Edge*>* elist = new list<Edge*>;
	for (it = leadingEdges_.begin(); it != leadingEdges_.end(); ++it) {
		Edge* edge = *it;
		for (int i = 0; i < 3; ++i) {
			Edge* twinedge = edge->getTwinEdge();
			// only one of the half-edges

			if ((twinedge == NULL && !skip_boundary_edges) ||
				(twinedge != NULL && ((size_t)edge >(size_t)twinedge)))
				elist->push_front(edge);

			edge = edge->getNextEdgeInFace();
		}
	}
	return elist;
}


bool IncrementalBlueNoise::InTriangle(Edge* edge, Node2d * target_node, int &mode)
{

//	Node2d* new_node = target_node;
	Node2d* n1 = edge->getSourceNode();
	Edge* e1 = edge;

	Edge* e2 = edge->getNextEdgeInFace();
	Node2d* n2 = e2->getSourceNode();

	Edge* e3 = e2->getNextEdgeInFace();
	Node2d* n3 = e3->getSourceNode();

	real ccw1 = CounterClockWise(n1, n2, target_node); // 如果target_node 在 1 2 的 右侧  则 顺时针， 负值
	real ccw2 = CounterClockWise(n2, n3, target_node); // 如果target_node 在 2 3 的 左侧  则 逆时针， 正值
	real ccw3 = CounterClockWise(n3, n1, target_node); 

	if (ccw1 > 0 && ccw2 > 0 && ccw3 > 0) 
	{
		mode = 0;
		return true;
	}
	else if (ccw1*ccw2*ccw3 == 0)
	{
		if (ccw2*ccw3 > 0)
			mode = 1;
		else if (ccw1*ccw3 > 0)
			mode = 2;
		else if (ccw1*ccw2 > 0)
			mode = 3;
		return true;
	}
	mode = -1;
	return false;
}

bool IncrementalBlueNoise::GlobalInsertNode(Node2d * target_node)
{


#ifdef _Debug_Demo
	flip_records.clear();
	split_records.clear();
#endif 

	//Locate the triangle
	list<Edge*>::const_iterator it;
	
	int mode = -1;

	for (it = leadingEdges_.begin(); it != leadingEdges_.end(); ++it)
	{
		Edge* edge = *it;
		bool  inTri = InTriangle(edge, target_node, mode);
		if (inTri)
		{
			if (mode == 0)
			{
				splitTriangle(*edge, *target_node);
				return true;
			}

			else if (mode == 1)
			{
				splitOnEdge(*edge, *target_node);
				return true;
			}

			else if (mode == 2)
			{
				edge = edge->getNextEdgeInFace();
				splitOnEdge(*edge, *target_node);
				return true;
			}
			
			else if (mode == 3)
			{
				edge = edge->getNextEdgeInFace();
				edge = edge->getNextEdgeInFace();
				splitOnEdge(*edge, *target_node);
				return true;
			}			

		}
	}
	return false;

}

/*          ________ vh_oppo
*          /\      /
*         /  \    /
*        /    \  /
*    vh /_<____\/
*/

Edge* IncrementalBlueNoise::splitOnEdge(Edge& edge, Node2d& point)
{


	Node2d* new_node = &point;
	Node2d* n1 = edge.getSourceNode();

	Edge* e1 = &edge;

	Edge* e2 = edge.getNextEdgeInFace();
	Node2d* n2 = e2->getSourceNode();

	Edge* e3 = e2->getNextEdgeInFace();
	Node2d* n3 = e3->getSourceNode();

	Edge* e1_n = new Edge;
	Edge* en_2 = new Edge;
	Edge* en_3 = new Edge;
	Edge* e3_n = new Edge;

	/*
	理顺关系
	*/
	e1_n->setSourceNode(n1);
	en_2->setSourceNode(new_node);
	en_3->setSourceNode(new_node);
	e3_n->setSourceNode(n3);

	en_3->setTwinEdge(e3_n);
	e3_n->setTwinEdge(en_3);

	e1_n->setNextEdgeInFace(en_3);
	en_3->setNextEdgeInFace(e3);
	e3->setNextEdgeInFace(e1_n);

	en_2->setNextEdgeInFace(e2);
	e2->setNextEdgeInFace(e3_n);
	e3_n->setNextEdgeInFace(en_2);

	Edge* leadingEdge;
	if (e1->isLeadingEdge())
		leadingEdge = e1;
	else if (e2->isLeadingEdge())
		leadingEdge = e2;
	else if (e3->isLeadingEdge())
		leadingEdge = e3;
	else
		return NULL;

	removeLeadingEdgeFromList(leadingEdge);

	addLeadingEdge(e1_n);
	addLeadingEdge(en_2);

	swapTestEdge(e2);
	swapTestEdge(e3);


	Edge* e1_twin = e1->getTwinEdge();

	if (e1_twin != NULL) {// if boundary .. 

		Node2d*  n1_twin = e1_twin->getSourceNode();

		Edge* e2_twin = e1_twin->getNextEdgeInFace();
		Node2d* n2_twin2 = e2_twin->getSourceNode();

		Edge* e3_twin = e2_twin->getNextEdgeInFace();
		Node2d* n3_twin = e3_twin->getSourceNode();



		Edge* e1_n_twin = new Edge;
		Edge* en_2_twin = new Edge;
		Edge* en_3_twin = new Edge;
		Edge* e3_n_twin = new Edge;
		/*
		理顺关系
		*/
		e1_n_twin->setSourceNode(n1_twin);
		en_2_twin->setSourceNode(new_node);
		en_3_twin->setSourceNode(new_node);
		e3_n_twin->setSourceNode(n3_twin);

		en_3_twin->setTwinEdge(e3_n_twin);
		e3_n_twin->setTwinEdge(en_3_twin);

		e1_n_twin->setNextEdgeInFace(en_3_twin);
		en_3_twin->setNextEdgeInFace(e3_twin);
		e3_twin->setNextEdgeInFace(e1_n_twin);

		en_2_twin->setNextEdgeInFace(e2_twin);
		e2_twin->setNextEdgeInFace(e3_n_twin);
		e3_n_twin->setNextEdgeInFace(en_2_twin);

		Edge* leadingEdge_twin;
		if (e1_twin->isLeadingEdge())
			leadingEdge_twin = e1_twin;
		else if (e2_twin->isLeadingEdge())
			leadingEdge_twin = e2_twin;
		else if (e3_twin->isLeadingEdge())
			leadingEdge_twin = e3_twin;
		else
			return NULL;

		removeLeadingEdgeFromList(leadingEdge_twin);

		/* special care */

		e1_n_twin->setTwinEdge(en_2);
		en_2->setTwinEdge(e1_n_twin);

		en_2_twin->setTwinEdge(e1_n);
		e1_n->setTwinEdge(en_2_twin);

		/* sepcial care end */

		addLeadingEdge(e1_n_twin);
		addLeadingEdge(en_2_twin);

		swapTestEdge(e2_twin);
		swapTestEdge(e3_twin);
	}

	return e1_n;
}
//--------------------------------------------------------------------------------------------------
Edge* IncrementalBlueNoise::splitTriangle(Edge& edge, Node2d& point) 
{

#ifdef _Debug_Demo
	split_records.push_back(point);
#endif 


	Node2d* new_node = &point;
	Node2d* n1 = edge.getSourceNode();
	Edge* e1 = &edge;

	Edge* e2 = edge.getNextEdgeInFace();
	Node2d* n2 = e2->getSourceNode();

	Edge* e3 = e2->getNextEdgeInFace();
	Node2d* n3 = e3->getSourceNode();


	Edge* e1_n = new Edge;
	Edge* e11_n = new Edge;
	Edge* e2_n = new Edge;
	Edge* e22_n = new Edge;
	Edge* e3_n = new Edge;
	Edge* e33_n = new Edge;

	e1_n->setSourceNode(n1);
	e11_n->setSourceNode(new_node);
	e2_n->setSourceNode(n2);
	e22_n->setSourceNode(new_node);
	e3_n->setSourceNode(n3);
	e33_n->setSourceNode(new_node);

	e1_n->setTwinEdge(e11_n);
	e11_n->setTwinEdge(e1_n);
	e2_n->setTwinEdge(e22_n);
	e22_n->setTwinEdge(e2_n);
	e3_n->setTwinEdge(e33_n);
	e33_n->setTwinEdge(e3_n);


	e1_n->setNextEdgeInFace(e33_n);
	e2_n->setNextEdgeInFace(e11_n);
	e3_n->setNextEdgeInFace(e22_n);

	e11_n->setNextEdgeInFace(e1);
	e22_n->setNextEdgeInFace(e2);
	e33_n->setNextEdgeInFace(e3);


	// and update old's next edge
	e1->setNextEdgeInFace(e2_n);
	e2->setNextEdgeInFace(e3_n);
	e3->setNextEdgeInFace(e1_n);



	// add the three new leading edges, 
	// Must remove the old leading edge from the list.
	// Use the field telling if an edge is a leading edge
	// NOTE: Must search in the list!!!


	Edge* leadingEdge;
	if (e1->isLeadingEdge())
		leadingEdge = e1;
	else if (e2->isLeadingEdge())
		leadingEdge = e2;
	else if (e3->isLeadingEdge())
		leadingEdge = e3;
	else
		return NULL;


	addLeadingEdge(e1_n);
	addLeadingEdge(e2_n);
	addLeadingEdge(e3_n);

	removeLeadingEdgeFromList(leadingEdge);


	
	swapTestEdge(e1);
	swapTestEdge(e2);
	swapTestEdge(e3);
	// Return a half edge incident to the new node (with the new node as source node)

	return e11_n;
}


//--------------------------------------------------------------------------------------------------
void IncrementalBlueNoise::swapTestEdge(Edge * diagonal) {

	// Note that diagonal is both input and output and it is always
	// kept in counterclockwise direction (this is not required by all 
	// finctions in ttl:: now)

	// Swap by rotating counterclockwise
	// Use the same objects - no deletion or new objects
	Edge* eR = diagonal;

	/*Code less*/
	if (eR == NULL) // Boundary 
		return;

	Edge* eL = eR->getTwinEdge();
	if (eL == NULL) // Boundary 
		return;
	/**/


	Edge* eL_1 = eL->getNextEdgeInFace();
	Edge* eL_2 = eL_1->getNextEdgeInFace();
	Edge* eR_1 = eR->getNextEdgeInFace();
	Edge* eR_2 = eR_1->getNextEdgeInFace();

	// avoid node to be dereferenced to zero and deleted
	Handle<Node2d> nR = eR_2->getSourceNode();
	Handle<Node2d> nL = eL_2->getSourceNode();

	//TEST HERE!
	/*Code less*/
	Handle<Node2d>  n1 = eL->getSourceNode();
	Handle<Node2d>  n2 = nR;
	Handle<Node2d>  n3 = eR->getSourceNode();
	Handle<Node2d>  n4 = nL;

	if (InCircle(n1, n2, n3, n4) >= 0) return;

#ifdef _Debug_Demo


	std::array<Handle<Node2d>, 4> flip_record;
	flip_record[0] = nR;
	flip_record[1] = n3;
	flip_record[2] = nL;
	flip_record[3] = n4;
	flip_records.push_back(flip_record);


#endif 

	
	/*Code less*/


	Edge* leL;
	if (eL->isLeadingEdge())
		leL = eL;
	else if (eL_1->isLeadingEdge())
		leL = eL_1;
	else if (eL_2->isLeadingEdge())
		leL = eL_2;


	Edge* leR;
	if (eR->isLeadingEdge())
		leR = eR;
	else if (eR_1->isLeadingEdge())
		leR = eR_1;
	else if (eR_2->isLeadingEdge())
		leR = eR_2;

	removeLeadingEdgeFromList(leL);
	removeLeadingEdgeFromList(leR);


	Edge* eUp = new Edge; // 
	Edge* eLow = new Edge; // 

	eUp->setSourceNode(nR.getPtr());

	eUp->setNextEdgeInFace(eL_2);
	eL_2->setNextEdgeInFace(eR_1);
	eR_1->setNextEdgeInFace(eUp);

	eLow->setSourceNode(nL.getPtr());

	eLow->setNextEdgeInFace(eR_2);
	eR_2->setNextEdgeInFace(eL_1);
	eL_1->setNextEdgeInFace(eLow);


	eUp->setTwinEdge(eLow);
	eLow->setTwinEdge(eUp);


	addLeadingEdge(eUp);
	addLeadingEdge(eLow);

	swapTestEdge(eL_2);
	swapTestEdge(eL_1);
}

