#pragma once
#include "calcubase.h"
#include "eigen/Eigen"
#include "OpenGLDC.h"
struct OBB
{
	vec_POINT3D theBox;//储存8个顶点
	//       4---7
	//      /|  /|
	//     5---6 |
	//     | 3-|-2
	//     |/  |/
	//     0---1
	//     */
	POINT3D theCenter;//中心点
	int ID;
	int preID;
	///////////////////////三个基向量
	VECTOR3D theX;
	VECTOR3D theY;
	VECTOR3D theZ;
	////////////////////////变换矩阵
	MATRIX3D MA;//世界坐标到局部坐标
	/////////////////////////	
	vec_PVERT* pBloop;
	vec_FACETTRI* pCUTface;
	BOOL isBranch;
	OBB()
	{
		theX.dx=0;theX.dy=0;theX.dz=0;
		theY.dx=0;theY.dy=0;theY.dz=0;
		theZ.dx=0;theZ.dy=0;theZ.dz=0;
		theCenter.x=0;theCenter.y=0;theCenter.z=0;
		ID = preID = -1;
		isBranch=FALSE;
		pBloop=NULL;
		pCUTface = NULL;
	}
};
class COBB
{
    public:
	static void makeOBB(vec_PFACETTRI &seg_BOX,OBB &ob);//创造OBB包围盒 
	static void updateBOX(vec_PFACETTRI &seg_BOX,OBB &ob);//更新包围盒
};
 