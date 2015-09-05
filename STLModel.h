#pragma once
#include "solidmesh.h"
#include "OpenGLDC.h"
#include "map";
#include "fstream";
#include <queue>
#include "eigen/Eigen"
#include "calcubase.h"
#include "solidPolygon.h"
#include "OBB.h"
#include <time.h>
using namespace Eigen;
using Eigen::MatrixXd;
using namespace std;

class CSTLModel :public CSolidMeshTri
{

public:
	CSTLModel(void);
	~CSTLModel(void);
	BOX3D m_one;
	long maxz0;
	long maxz1;
	long mazx0;
	long maxx1;
	Ccalcubase* pCalcubase;
public:
	BOOL LoadSTLFile(LPCTSTR stlfile);            //读取STL文件
	virtual void Draw(COpenGLDC* pDC);             //显示函数
	virtual void Draw1(COpenGLDC* pDC,vector<vector<PFACETTRI>>& Rbox);
	virtual void Draw3(COpenGLDC* pDC,vector<vector<PFACETTRI>>& Rbox);
	virtual void Draw2(COpenGLDC* pDC);
	virtual void Draw_obb(COpenGLDC* pDC,vector<OBB>& vec_obb);  
	virtual void Draw_ccPath(COpenGLDC* pDC,vec_PPOLYPOLYGON&	m_vecPPolyPolygons); 
	virtual bool	IsEmpty()	;
	virtual	void	UpdateBox()	;
	virtual void    GetBox(double x0,double y0,double z0,double x1,double y1,double z1);
	vector<POINT3D>  vec_box; 
};