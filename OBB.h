#pragma once
#include "calcubase.h"
#include "eigen/Eigen"
#include "OpenGLDC.h"
struct OBB
{
	vec_POINT3D theBox;//����8������
	//       4---7
	//      /|  /|
	//     5---6 |
	//     | 3-|-2
	//     |/  |/
	//     0---1
	//     */
	POINT3D theCenter;//���ĵ�
	int ID;
	int preID;
	///////////////////////����������
	VECTOR3D theX;
	VECTOR3D theY;
	VECTOR3D theZ;
	////////////////////////�任����
	MATRIX3D MA;//�������굽�ֲ�����
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
	static void makeOBB(vec_PFACETTRI &seg_BOX,OBB &ob);//����OBB��Χ�� 
	static void updateBOX(vec_PFACETTRI &seg_BOX,OBB &ob);//���°�Χ��
};
 