#pragma once 
#include "calcubase.h"
#include "OpenGLDC.h"
class CSDF
{
public:

	/////////////////////////////////////////////////////SDF算法
	/////////////////////////////////////////////////SDF聚类
	double SDF_calcu(PFACETTRI faceA,Octree* pOctree,vec_PFACETTRI& m_vecPFacetTri);//计算一个面片的SDF值
	void SDF_kmeans(int k,vector<vector<PVERT>> &RBLOOP,vec_PFACETTRI& m_vecPFacetTri,vector<POINT3D>& vec_box,vector<vector<PFACETTRI>>& Rbox);//k_means聚类
	void SDF_gram(int k,vector<vector<PVERT>> &RBLOOP,vector<double>& sdf_gram,vec_PFACETTRI& m_vecPFacetTri,vector<POINT3D>& vec_box,vec_PVERT& m_vecPVert);//SDF直方图
	void SDF_GMM(int k,vector<vector<PVERT>> &RBLOOP,vec_PFACETTRI& m_vecPFacetTri,vector<POINT3D>& vec_box,vector<vector<PFACETTRI>>& Rbox);//GMM聚类
	///////////////////////////////////////////////////////新特征
	double SDF_featurPOINT(PVERT& Apot,Octree* pOctree,vec_PVERT& sdf_point,vec_PFACETTRI& m_vecPFacetTri);
	double SDF_featurFACE(PFACETTRI faceA,Octree* pOctree,vec_PFACETTRI& m_vecPFacetTri);
	void SDF_getPOINT(vec_PFACETTRI& m_vecPFacetTri,vec_PVERT& m_vecPVert,vector<POINT3D>& vec_box);
	//vector<vector<PFACETTRI>> Rbox;
protected:
private:
};