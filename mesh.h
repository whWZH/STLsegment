#pragma once

#include "GeomBase.h"


struct   HalfEdge;
struct   Edge;
struct   Vertex;
struct	 Facet;

typedef struct   HalfEdge
{
	Vertex*   pVertEnd;   // vertex at the end of the half-edge
	HalfEdge* pHEdgePair;   // oppositely oriented adjacent half-edge 
	Facet*    pFacetAdj;         // face the half-edge borders
	HalfEdge* pHEdgeNext;   // next half-edge around the face
	HalfEdge* pHEdgePre;
	BOOL      bStatus;
	BOOL      bused;
	BOOL      be_pair;
	HalfEdge(){
		pVertEnd = 0;
		pHEdgePair = 0;
		pFacetAdj = 0;
		pHEdgeNext = 0;
		pHEdgePre = 0;
		bused = 0;
		be_pair=0;
		bStatus=FALSE;
	};
}HEDGE,*PHEDGE;


typedef struct   Edge
{
	Vertex*   pVertStart;
	Vertex*   pVertEnd;
	HalfEdge* pHEdgeAdj;
	BOOL      bStatus;

	Edge(){
		pVertStart = 0;
		pVertEnd = 0;
		pHEdgeAdj = 0;
		bStatus=FALSE;
	};
}EDGE, *PEDGE;


typedef struct   Vertex: public tagPoint3D
{
	HalfEdge* pHEdgeOut;  // one of the half-edges emanating from the vertex
	BOOL      bStatus;
	BOOL      bused;
	BOOL      bSharpVer;
	BOOL      theEND;
	int       ID;
	VECTOR3D  Normal;
	Vertex(){pHEdgeOut = 0;bStatus=FALSE;bused=0;bSharpVer = FALSE;theEND=0;ID=0;Normal.dx=Normal.dy=Normal.dz=0;};
}VERT, *PVERT;

typedef struct	Facet
{
	HalfEdge* pHEdge;  // one of the half-edges bordering the face
	BOOL      bStatus;
	Facet(){pHEdge = 0;bStatus=FALSE; };
} FACET, *PFACET;

typedef struct FacetTri:public Facet
{
	PVERT     m_PVerts[3];
	PVECTOR3D m_PFacetNorm;
	bool      becut;
	int       ID;
	FacetTri()
	{
		m_PVerts[0] = m_PVerts[1] =m_PVerts[2] = 0;
		m_PFacetNorm =0;
		becut=0;
		ID=-1;
	};
}FACETTRI, *PFACETTRI;


typedef	vector<PHEDGE> vec_PHEDGE;
typedef vec_PHEDGE::iterator it_vec_PHEDGE;
typedef	vector<HEDGE> vec_HEDGE;
typedef vec_HEDGE::iterator it_vec_HEDGE;

typedef	vector<PEDGE> vec_PEDGE;
typedef vec_PEDGE::iterator it_vec_PEDGE;
typedef	vector<EDGE> vec_EDGE;
typedef vec_EDGE::iterator it_vec_EDGE;

typedef	vector<PVERT> vec_PVERT;
typedef vec_PVERT::iterator it_vec_PVERT;
typedef	vector<VERT> vec_VERT;
typedef vec_VERT::iterator it_vec_VERT;

typedef	vector<PFACET> vec_PFACET;
typedef vec_PFACET::iterator it_vec_PFACET;
typedef	vector<FACET> vec_FACET;
typedef vec_FACET::iterator it_vec_FACET;

typedef	vector<PFACETTRI> vec_PFACETTRI;
typedef vec_PFACETTRI::iterator it_vec_PFACETTRI;
typedef	vector<FACETTRI> vec_FACETTRI;
typedef vec_FACETTRI::iterator it_vec_FACETTRI;
