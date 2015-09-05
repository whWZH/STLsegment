#if !defined (OCTREE_H__INCLUDE_)
#define OCTREE_H__INCLUDE_
#include "calcubase.h"
typedef double* Vec3;
extern Vec3 makeVec3( double x, double y, double z);

/* ------------------------------------------------------------------------ */
class pvex_nor
{
public:
	pvex_nor(double xx,double yy,double zz,int id){x=xx;y=yy;z=zz;ID=id;};
	double x,y,z;
	int ID;
};
typedef struct {
   Vec3 p[3];
} TRIANGLE;

typedef struct {
	Vec3 p[8];
	double val[8];
} GRIDCELL;
class Octree
{
public:
	Vec3 min;  /* bounds */
	Vec3 max;
	double value[8];
	int density;
	char at_max_depth;
	char not_fully_divided;
	CPtrList vex;
	CPtrList normal;
	Octree** children;
    Octree* make_octree( Vec3 min, Vec3 max );
    void isoface(Octree* o);
	void split_octree( Octree* o );
    int subdivide_octree( int min_depth, int max_depth, Octree* o );
	void subpoint( Octree* o,int oc,Vec3 min, Vec3 max);
    int octree_needs_to_be_split( Octree* o );
    void marchingcube(int depth,Octree *o);
	int recursively_evaluate_octree( int min_depth, int max_depth, int current_depth, Octree* o);
    double evaluate_point( Vec3 pos,Octree *o );
    int evaluate1_point( Octree* o);
    int octree_needs_to_be_subdivided( Octree* o );
	///////////////////////////////////////////////////
	///////////////////////////////////////////////八叉树相交判断
	void creat_octree(Octree*& pOctree,vector<PFACETTRI> &vec_face,vector<POINT3D>&  vec_box);//建立八叉树
	void delet_octree(Octree*& pOctree);//析构八叉树
	void delet_octree_1(Octree*& pOctree);
	PFACETTRI get_face_new(Octree* &pOctree,VECTOR3D LineNor,POINT3D orig,PFACETTRI faceA,vec_PFACETTRI& m_vecPFacetTri);//回溯法遍历
	PFACETTRI get_face_new1(Octree* &pOctree,VECTOR3D LineNor,POINT3D orig,PFACETTRI faceA,vec_PFACETTRI& m_vecPFacetTri);//回溯法遍历
	BOOL get_face(Octree* pOctree,VECTOR3D LineNor,POINT3D orig,PFACETTRI faceA,vector<int>& ID_crossface);//求与射线相交的三角面片,树、面片集、射线、box
	BOOL cros_BigCube(Octree* pOctree,VECTOR3D LineNor,POINT3D orig);//判断射线是否跟父节点所在的大三角相交
	
protected:
private:
};
//typedef struct _Octree
//{
//  Vec3 min;  /* bounds */
//  Vec3 max;
//
//  double value[8];
//  int density;
//  char at_max_depth;
//  char not_fully_divided;
//  CPtrList vex;
//  CPtrList normal;
//  struct _Octree** children;
//
//} Octree;
//extern Octree* make_octree( Vec3 min, Vec3 max );
//extern void isoface(Octree* o);
//extern int subdivide_octree( int min_depth, int max_depth, Octree* o );
//extern int octree_needs_to_be_split( Octree* o );
//void marchingcube(int depth,Octree *o);
//
///* the first function of interest is: 
//
//   evaluatePoint( ) 
//   
//   which takes a 3d point and returns a double value indicating which feature 
//   this point is nearet to.
//
//   see the custom.c file for some examples
//
//*/
//
//extern double evaluate_point( Vec3 pos,Octree *o );
//extern int evaluate1_point( Octree* o);
///* the other function of interest is: 
//
//   octreeNeedToBeSplit( ) 
//   
//   which determines whether subdivision needs to take place.
//   
//   it should examine the 8 corner values and
//   return 1 to subdivide and 0 to not subdivide.
//
//   see the custom.c file for some examples
//
//*/
//
//extern int octree_needs_to_be_subdivided( Octree* o );
//
//
///* ------------------------------------------------------------------------ */

#endif
