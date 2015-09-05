#pragma once 
#include "calcubase.h"
#include "OpenGLDC.h"
class CSDF
{
public:

	/////////////////////////////////////////////////////SDF�㷨
	/////////////////////////////////////////////////SDF����
	double SDF_calcu(PFACETTRI faceA,Octree* pOctree,vec_PFACETTRI& m_vecPFacetTri);//����һ����Ƭ��SDFֵ
	void SDF_kmeans(int k,vector<vector<PVERT>> &RBLOOP,vec_PFACETTRI& m_vecPFacetTri,vector<POINT3D>& vec_box,vector<vector<PFACETTRI>>& Rbox);//k_means����
	void SDF_gram(int k,vector<vector<PVERT>> &RBLOOP,vector<double>& sdf_gram,vec_PFACETTRI& m_vecPFacetTri,vector<POINT3D>& vec_box,vec_PVERT& m_vecPVert);//SDFֱ��ͼ
	void SDF_GMM(int k,vector<vector<PVERT>> &RBLOOP,vec_PFACETTRI& m_vecPFacetTri,vector<POINT3D>& vec_box,vector<vector<PFACETTRI>>& Rbox);//GMM����
	///////////////////////////////////////////////////////������
	double SDF_featurPOINT(PVERT& Apot,Octree* pOctree,vec_PVERT& sdf_point,vec_PFACETTRI& m_vecPFacetTri);
	double SDF_featurFACE(PFACETTRI faceA,Octree* pOctree,vec_PFACETTRI& m_vecPFacetTri);
	void SDF_getPOINT(vec_PFACETTRI& m_vecPFacetTri,vec_PVERT& m_vecPVert,vector<POINT3D>& vec_box);
	//vector<vector<PFACETTRI>> Rbox;
protected:
private:
};