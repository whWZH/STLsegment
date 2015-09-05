#include "stdafx.h"
#include "calcubase.h"
////////////////////////////驱动平面函数////////////////////////////////
CPlane::CPlane(void)
{
	m_PointIn = POINT3D(0,0,0);
	m_Normal  = VECTOR3D(0,0,1);
	a = 0;b = 0;c = 1;d = 0;
}

CPlane::CPlane(POINT3D point,VECTOR3D normal)
{
	m_PointIn = point;
	m_Normal  = normal.GetNormal();
	a = m_Normal.dx; b = m_Normal.dy; c = m_Normal.dz;
	d = m_Normal|(POINT3D(0,0,0)-m_PointIn);             //等于原点到平面的距离
}

CPlane::~CPlane(void)
{
}

void CPlane::GetParameter()
{
	m_Normal.Normalize();
	a = m_Normal.dx; b = m_Normal.dy; c = m_Normal.dz;
	d = m_Normal|(POINT3D(0,0,0)-m_PointIn);             //等于原点到平面的距离
}
// 计算线段与平面的交点
BOOL Ccalcubase::IntersectLinePlane(PPOINT3D pStartP, PPOINT3D pEndP, CPlane Plane, POINT3D& InterP)
{
	VECTOR3D vSE = *pEndP - *pStartP;
	Plane.GetParameter();
	//Plane.d   = -((Plane.m_PointIn-POINT3D(0,0,0))|Plane.m_Normal);
	double ts = ((*pStartP-POINT3D(0,0,0))|Plane.m_Normal)+Plane.d;
	double te = ((*pEndP-POINT3D(0,0,0))|Plane.m_Normal)+Plane.d;
	double t = fabs(ts)+fabs(te);

	if (ts*te<0||IS_ZERO(ts*te))                    //两点在面的两侧，或在面上
	{
		if (IS_ZERO(ts*te))                //有点在面上
		{
			if (IS_ZERO(ts))
			{
				InterP = *pStartP;    //第一个点在面上或两个点都在面上，第一个点为交点
				return TRUE;
			} 
			else
			{
				InterP = *pEndP;     //第二个点在面上，第二个点为交点
				return TRUE;
			}
		} 
		else                        //点在面的两侧
		{
			InterP = *pStartP + vSE*fabs(ts)/t;
			return TRUE;
		}
	} 
	else
	{
		return FALSE;         //无交点
	}
	//return 0;
}
BOOL Ccalcubase::IntersectLinePlane(PHEDGE pHEdge, CPlane Plane, POINT3D& InterP)
{
	PPOINT3D pStartP, pEndP;
	pStartP = pHEdge->pHEdgePre->pVertEnd;
	pEndP   = pHEdge->pVertEnd;
	VECTOR3D vSE = *pEndP - *pStartP;
	Plane.GetParameter();
	//Plane.d   = -((Plane.m_PointIn-POINT3D(0,0,0))|Plane.m_Normal);
	double ts = ((*pStartP-POINT3D(0,0,0))|Plane.m_Normal)+Plane.d;
	double te = ((*pEndP-POINT3D(0,0,0))|Plane.m_Normal)+Plane.d;
	double t = fabs(ts)+fabs(te);

	if (ts*te<0||IS_ZERO(ts*te))                    //两点在面的两侧，或在面上
	{
		if (IS_ZERO(ts*te))                //有点在面上
		{
			if (IS_ZERO(te))
			{
				InterP = *pEndP;     //第二个点在面上或两个点都在面上，第二个点为交点				
				return TRUE;
			} 
			else
			{
				InterP = *pStartP;    //第一个点在面上，第一个点为交点
				return TRUE;
			}
			//return FALSE;
		} 
		else                        //点在面的两侧
		{
			InterP = *pStartP + vSE*fabs(ts)/t;
			return TRUE;
		}
	} 
	else
	{
		return FALSE;         //无交点
	}
	//return 0;
}
BOOL Ccalcubase::IntersectLinePlane(POINT3D LinePt,VECTOR3D LineNor, PFACETTRI pFacTri, POINT3D& InterP)
{
	double t;             //直线参数方程x = LinePt.x+ LineNor.x * t; y = LinePt.y+ LineNor.y * t; z = LinePt.z+ LineNor.z * t
	LineNor.Normalize();
	POINT3D Pt;            //点法式平面
	VECTOR3D FacNor;
	Pt = (POINT3D)*(pFacTri->m_PVerts[0]);
	FacNor = *(pFacTri->m_PFacetNorm);

	POINT3D InterP_temp;  //平面交点
	VECTOR3D V01,V12,V20,V0i,V1i,V2i;
	BOOL  IsInFac = FALSE;    //判断是否在三角形内

	if (!IS_ZERO(LineNor|FacNor))
	{
		t = ((Pt-LinePt)|FacNor)/(LineNor|FacNor);   //t = ((Pt.x – LinePt.x)*FacNor.x+(Pt.y – LinePt.y)*FacNor.y+(Pt.z – LinePt.z)*FacNor.z) 
		//   /(FacNor.x* LineNor.x+ FacNor.y* LineNor.y+ FacNor.z* LineNor.z)    
		InterP_temp = LinePt + LineNor*t;       //将求得参数t带入直线参数方程求交点

		V01 = (POINT3D)*(pFacTri->m_PVerts[1]) - (POINT3D)*(pFacTri->m_PVerts[0]);
		V12 = (POINT3D)*(pFacTri->m_PVerts[2]) - (POINT3D)*(pFacTri->m_PVerts[1]);
		V20 = (POINT3D)*(pFacTri->m_PVerts[0]) - (POINT3D)*(pFacTri->m_PVerts[2]);
		V0i = InterP_temp - (POINT3D)*(pFacTri->m_PVerts[0]);
		V1i = InterP_temp - (POINT3D)*(pFacTri->m_PVerts[1]);
		V2i = InterP_temp - (POINT3D)*(pFacTri->m_PVerts[2]);

		if (((V01*V0i).dz>0||IS_ZERO((V01*V0i).dz))&&((V12*V1i).dz>0||IS_ZERO((V12*V1i).dz))&&((V20*V2i).dz>0||IS_ZERO((V20*V2i).dz))) IsInFac = TRUE;    //若叉乘的符号相同则点在面片上
		if (((V01*V0i).dz<0||IS_ZERO((V01*V0i).dz))&&((V12*V1i).dz<0||IS_ZERO((V12*V1i).dz))&&((V20*V2i).dz<0||IS_ZERO((V20*V2i).dz))) IsInFac = TRUE;

		if (IsInFac)
		{
			InterP = InterP_temp;  //平面上交点为面片上交点
			return TRUE;
		}
	}
	return FALSE;         //无交点
}
//三角面片与平面交点，2013.7.18
//其中pHEdge是前一个三角面片半边，InterpHE是当前正在求交的三角面片的半边
BOOL Ccalcubase::IntersectFacPlane(POINT3D PointPre,PHEDGE pHEdge,POINT3D Point, CPlane Plane,PHEDGE& InterpHE, POINT3D& InterP)
{                                //（前前一交点，前一交点所在半边，前一交点，截交面，所求交点所在半边，所求交点）
	//第一个参数是为了判断后一个交点是网格顶点的情况
	PHEDGE pHE_Inter;  //循环终止的半边
	BOOL   bInter;     //判断是否相交

	if (*(pHEdge->pVertEnd) ==Point||*(pHEdge->pHEdgePair->pVertEnd) ==Point)//判断交点是否为网格顶点	
	{
		//当上一个交点是网格顶点时
		if (*(pHEdge->pVertEnd) ==Point) InterpHE = pHEdge->pHEdgePair->pHEdgeNext;
		else                             InterpHE = pHEdge->pHEdgeNext;

		pHE_Inter = InterpHE;
		do      //从一阶领域外围半边判断交点
		{
			bInter = Ccalcubase::IntersectLinePlane(InterpHE,Plane,InterP);
			if (bInter == TRUE)
			{
				if (InterP == PointPre) bInter  = FALSE; //判断交点是否为前交点
				else break;
			}
			InterpHE = InterpHE->pHEdgeNext->pHEdgePair->pHEdgeNext;
		}while(InterpHE != pHE_Inter);

	}
	else
	{
		InterpHE = pHEdge->pHEdgePair->pHEdgeNext;
		bInter = Ccalcubase::IntersectLinePlane(InterpHE,Plane,InterP);
		if (bInter == FALSE)
		{
			InterpHE = InterpHE->pHEdgeNext;
			bInter = Ccalcubase::IntersectLinePlane(InterpHE,Plane,InterP);
		}
	}
	return bInter;
}
//求空间点是否在面片上,2014.1.4
BOOL Ccalcubase::IsPointInFac(POINT3D Pt,PFACETTRI pFacTri)
{
	double xx0 = pFacTri->m_PVerts[0]->x,
		yy0 = pFacTri->m_PVerts[0]->y,
		xx1 = pFacTri->m_PVerts[0]->x,
		yy1 = pFacTri->m_PVerts[0]->y,
		zz0 = pFacTri->m_PVerts[0]->z,
		zz1 = pFacTri->m_PVerts[0]->z;

	for (int i=1;i<3;i++)
	{
		xx0 = (xx0<pFacTri->m_PVerts[i]->x)?xx0:pFacTri->m_PVerts[i]->x;
		xx1 = (xx1>pFacTri->m_PVerts[i]->x)?xx1:pFacTri->m_PVerts[i]->x;
		yy0 = (yy0<pFacTri->m_PVerts[i]->y)?yy0:pFacTri->m_PVerts[i]->y;
		yy1 = (yy1>pFacTri->m_PVerts[i]->y)?yy1:pFacTri->m_PVerts[i]->y;
		zz0 = (zz0<pFacTri->m_PVerts[i]->z)?zz0:pFacTri->m_PVerts[i]->z;
		zz1 = (zz1>pFacTri->m_PVerts[i]->z)?zz1:pFacTri->m_PVerts[i]->z;
	}

	if (Pt.x>=xx0 && Pt.x<=xx1 && Pt.y>=yy0 && Pt.y<=yy1 && Pt.z>=zz0 && Pt.z<=zz1)
	{
		VECTOR3D V0,V1,V2;
		double   A0,A1,A2;
		V0 = *(pFacTri->m_PVerts[0]) - Pt;
		V1 = *(pFacTri->m_PVerts[1]) - Pt;
		V2 = *(pFacTri->m_PVerts[2]) - Pt;
		A0 = _AngleBetween3D(V0,V2);
		A1 = _AngleBetween3D(V1,V2);
		A2 = _AngleBetween3D(V1,V0);

		if (IS_ZERO(2*PI-A0-A1-A2))	return TRUE;
		else                        return FALSE;
	} 
	else return FALSE;

}
/////////////////////////////////建立局部坐标系
Ccalcubase::Ccalcubase(void)
{

}
Ccalcubase::~Ccalcubase(void)
{

}
//整个模型上,两个点之间的最短路径
void Ccalcubase::Dijkstra_Point_ST_EN_new(PVERT ST,PVERT EN,map<PVERT,PVERT>& Path,vector<PVERT>& BLOOP,vec_PVERT& m_vecPVert)
{
	queue<PVERT> Q; 
	map<PVERT,double> Dist; 
	vector<PVERT> Pnear;
	for (int i=0;i<m_vecPVert.size();i++)
	{
		m_vecPVert[i]->bStatus=0;
		Dist.insert(pair<PVERT,double>(m_vecPVert[i],999));
	}
	Dist[ST]=0; 
	Q.push(ST);
	Path[ST]=ST->pHEdgeOut->pVertEnd;
	ST->bStatus=1;
	while(!Q.empty())
	{
		vec_PVERT vecpVer;
		FindOneRing(Q.front(),vecpVer);
		for (int i=0;i<vecpVer.size();i++)
		{
			Pnear.push_back(vecpVer[i]);
			double k=1;//
			if ((Dist[Q.front()]+k)<Dist[vecpVer[i]])
			{
				Path[vecpVer[i]]=Q.front();
				Dist[vecpVer[i]]=Dist[Q.front()]+k;
			}
		}
		while(!Pnear.empty())
		{
			if (Pnear.size()>1)
			{
				for (int i=0;i<Pnear.size()-1;i++)
				{
					if (Dist[Pnear[i]]<Dist[Pnear[i+1]])
					{
						swap(Pnear[i],Pnear[i+1]);
					}
				}
				if (Pnear.back()->bStatus==0)
				{
					Pnear.back()->bStatus=1;
					Q.push(Pnear.back());
				}
				Pnear.pop_back();
			}
			if (Pnear.size()==1)
			{
				if (Pnear.back()->bStatus==0)
				{
					Pnear.back()->bStatus=1;
					Q.push(Pnear.back());
				}
				Pnear.pop_back();
			}
		}
		Q.pop();
		if (Q.size()>0)
		{
			if (Q.front()==EN)
			{
				break;
			}
		}
	}
}
///////////////////////////////射线求交点
BOOL Ccalcubase::IntersectTriangle(POINT3D orig,VECTOR3D dir,POINT3D v0,POINT3D v1,POINT3D v2,POINT3D &JPOT)
{
	double t,u,v;
	VECTOR3D E1 = v1 - v0;
	VECTOR3D E2 = v2 - v0;
	VECTOR3D P = dir*E2; 
	double det = E1|P;
	VECTOR3D T;
	VECTOR3D QJ;
	if( det >0 )
	{
		T = orig - v0;
	}
	else
	{
		T = v0 - orig;
		det = -det;
	}

	if (det < 0.0001)
		return false;
	u = T|P;
	if( u < 0 || u > det )
		return false;
	QJ=(T*E1);
	v = (dir|QJ);
	if( v < 0 || u + v > det )
		return false;
	t = (E2|QJ);
	if (t<0)
	{
		return false;
	}
	double fInvDet = (double)1 /(double)det;
	t = fInvDet;
	u = fInvDet;
	v  = fInvDet;
	JPOT=v0*(1 - u - v) + v1*u + v2*v;
	return true;
}
BOOL Ccalcubase::IntersectTriangle(POINT3D orig,VECTOR3D dir,POINT3D v0,POINT3D v1,POINT3D v2)
{
	double t,u,v;
	VECTOR3D E1 = v1 - v0;
	VECTOR3D E2 = v2 - v0;
	VECTOR3D P = dir*E2; 
	double det = E1|P;
	VECTOR3D T;
	VECTOR3D QJ;
	if( det >0 )
	{
		T = orig - v0;
	}
	else
	{
		T = v0 - orig;
		det = -det;
	}

	if( det < 0)
		return false;
	u = T|P;
	if( u < 0 || u > det )
		return false;
	QJ=(T*E1);
	v = (dir|QJ);
	if( v < 0 || u + v > det )
		return false;
	t = (E2|QJ);
	if (t<0)
	{
		return false;
	}
	double fInvDet = 1 / det;
	t = fInvDet;
	u = fInvDet;
	v  = fInvDet;
	return TRUE;
}
void Ccalcubase::create_current_Frenet(MATRIX3D& currentM,VECTOR3D LineNor,POINT3D Center)
{
	double cosT,cosU,cosV;
	cosT=cosU=cosV=0;
	VECTOR3D X,Y,Z;
	double A,B,C,D;
	A=B=C=D=0;
	X.dx=1;X.dy=0;X.dz=0;Y.dx=0;Y.dy=1;Y.dz=0;Z.dx=0;Z.dy=0;Z.dz=1;
	cosT=(LineNor|X);
	cosU=(LineNor|Y);
	cosV=(LineNor|Z);
	D=sqrt(cosT*cosT+cosU*cosU);
	A=(-Center.x*cosT*cosV/D)-(Center.y*cosU*cosV/D)+Center.z*D;
	B=(Center.x*cosU/D)-(Center.y*cosT/D);
	C=-Center.x*cosT-Center.y*cosU-Center.z*cosV;
	currentM.A[0][0]=cosT*cosV/D;currentM.A[0][1]=(-cosU/D);currentM.A[0][2]=cosT;currentM.A[0][3]=0;
	currentM.A[1][0]=cosU*cosV/D;currentM.A[1][1]=cosT/D;currentM.A[1][2]=cosU;currentM.A[1][3]=0;
	currentM.A[2][0]=(-D);currentM.A[2][1]=0;currentM.A[2][2]=cosV;currentM.A[2][3]=0;
	currentM.A[3][0]=A;currentM.A[3][1]=B;currentM.A[3][2]=C;currentM.A[3][3]=1;
}
void Ccalcubase::create_BACK_Frenet(MATRIX3D& currentB,VECTOR3D LineNor,POINT3D Center)
{
	double cosT,cosU,cosV;
	cosT=cosU=cosV=0;
	VECTOR3D X,Y,Z;
	X.dx=1;X.dy=0;X.dz=0;Y.dx=0;Y.dy=1;Y.dz=0;Z.dx=0;Z.dy=0;Z.dz=1;
	cosT=(LineNor|X);
	cosU=(LineNor|Y);
	cosV=(LineNor|Z);
	double ABC=sqrt(cosT*cosT+cosU*cosU);
	currentB.A[0][0]=(cosT*cosV)/ABC;currentB.A[0][1]=(cosU*cosV)/ABC;currentB.A[0][2]=-ABC;currentB.A[0][3]=0;
	currentB.A[1][0]=(-cosU)/ABC;currentB.A[1][1]=cosT/ABC;currentB.A[1][2]=0;currentB.A[1][3]=0;
	currentB.A[2][0]=cosT;currentB.A[2][1]=cosU;currentB.A[2][2]=cosV;currentB.A[2][3]=0;
	currentB.A[3][0]=Center.x;currentB.A[3][1]=Center.y;currentB.A[3][2]=Center.z;currentB.A[3][3]=1;
}
///////////////////////////////建立绕X轴旋转的坐标
void Ccalcubase::create_current_matriX(MATRIX3D& matriX,double a)
{
	/*MATRIX3D matriX;
	double m[4];*/
	a=(a/180)*3.1415926;
	matriX.A[0][0]=cos(a);matriX.A[0][1]=0;matriX.A[0][2]=(-sin(a));matriX.A[0][3]=0;
	matriX.A[1][0]=0;matriX.A[1][1]=1;matriX.A[1][2]=0;matriX.A[1][3]=0;
	matriX.A[2][0]=sin(a);matriX.A[2][1]=0;matriX.A[2][2]=cos(a);matriX.A[2][3]=0;
	matriX.A[3][0]=0;matriX.A[3][1]=0;matriX.A[3][2]=0;matriX.A[3][3]=1;
	//m[0]=LineNor.dx;m[1]=LineNor.dy;m[2]=LineNor.dz;
	//m[3]=1;
	//LineNor.dx=(matriX.A[0][0]*m[0]+matriX.A[1][0]*m[1]+matriX.A[2][0]*m[2]+matriX.A[3][0]*m[3]);
	//LineNor.dy=(matriX.A[0][1]*m[0]+matriX.A[1][1]*m[1]+matriX.A[2][1]*m[2]+matriX.A[3][1]*m[3]);
	//LineNor.dz=(matriX.A[0][2]*m[0]+matriX.A[1][2]*m[1]+matriX.A[2][2]*m[2]+matriX.A[3][2]*m[3]);
}
void Ccalcubase::create_current_matriZ(MATRIX3D& matriZ,double a)
{
	a=(a/180)*3.1415926;
	matriZ.A[0][0]=cos(a);matriZ.A[0][1]=sin(a);matriZ.A[0][2]=0;matriZ.A[0][3]=0;
	matriZ.A[1][0]=(-sin(a));matriZ.A[1][1]=cos(a);matriZ.A[1][2]=0;matriZ.A[1][3]=0;
	matriZ.A[2][0]=0;matriZ.A[2][1]=0;matriZ.A[2][2]=1;matriZ.A[2][3]=0;
	matriZ.A[3][0]=0;matriZ.A[3][1]=0;matriZ.A[3][2]=0;matriZ.A[3][3]=1;
}
void Ccalcubase::create_cone(vector<VECTOR3D>& m_vecVECTOR3D,VECTOR3D LineNor,POINT3D orig,double angA,double angB)
{
	m_vecVECTOR3D.clear();
	VECTOR3D tempLine;
	MATRIX3D matriX,matriZ,currentB,currentM;
	double m[4];
	create_current_Frenet(currentM,LineNor,orig);
	create_BACK_Frenet(currentB,LineNor,orig);
	create_current_matriX(matriX,angA/3);
	create_current_matriZ(matriZ,angB);
	m[0]=LineNor.dx;m[1]=LineNor.dy;m[2]=LineNor.dz;
	m[3]=1;
	LineNor.dx=(currentM.A[0][0]*m[0]+currentM.A[1][0]*m[1]+currentM.A[2][0]*m[2]+currentM.A[3][0]*m[3]);
	LineNor.dy=(currentM.A[0][1]*m[0]+currentM.A[1][1]*m[1]+currentM.A[2][1]*m[2]+currentM.A[3][1]*m[3]);
	LineNor.dz=(currentM.A[0][2]*m[0]+currentM.A[1][2]*m[1]+currentM.A[2][2]*m[2]+currentM.A[3][2]*m[3]);
	tempLine.dx=0;tempLine.dy=0;tempLine.dz=1;
	m[0]=tempLine.dx;m[1]=tempLine.dy;m[2]=tempLine.dz;
	m[3]=1;
	tempLine.dx=(matriX.A[0][0]*m[0]+matriX.A[1][0]*m[1]+matriX.A[2][0]*m[2]+matriX.A[3][0]*m[3]);
	tempLine.dy=(matriX.A[0][1]*m[0]+matriX.A[1][1]*m[1]+matriX.A[2][1]*m[2]+matriX.A[3][1]*m[3]);
	tempLine.dz=(matriX.A[0][2]*m[0]+matriX.A[1][2]*m[1]+matriX.A[2][2]*m[2]+matriX.A[3][2]*m[3]);
	for (int i=0;i<(360/angB);i++)
	{
		double mm[4];
		mm[0]=tempLine.dx;mm[1]=tempLine.dy;mm[2]=tempLine.dz;
		mm[3]=1;
		tempLine.dx=(matriZ.A[0][0]*mm[0]+matriZ.A[1][0]*mm[1]+matriZ.A[2][0]*mm[2]+matriZ.A[3][0]*mm[3]);
		tempLine.dy=(matriZ.A[0][1]*mm[0]+matriZ.A[1][1]*mm[1]+matriZ.A[2][1]*mm[2]+matriZ.A[3][1]*mm[3]);
		tempLine.dz=(matriZ.A[0][2]*mm[0]+matriZ.A[1][2]*mm[1]+matriZ.A[2][2]*mm[2]+matriZ.A[3][2]*mm[3]);
		m_vecVECTOR3D.push_back(tempLine);
	}
	m[0]=tempLine.dx;m[1]=tempLine.dy;m[2]=tempLine.dz;
	m[3]=1;
	tempLine.dx=(matriX.A[0][0]*m[0]+matriX.A[1][0]*m[1]+matriX.A[2][0]*m[2]+matriX.A[3][0]*m[3]);
	tempLine.dy=(matriX.A[0][1]*m[0]+matriX.A[1][1]*m[1]+matriX.A[2][1]*m[2]+matriX.A[3][1]*m[3]);
	tempLine.dz=(matriX.A[0][2]*m[0]+matriX.A[1][2]*m[1]+matriX.A[2][2]*m[2]+matriX.A[3][2]*m[3]);
	for (int i=0;i<(360/angB);i++)
	{
		double mm[4];
		mm[0]=tempLine.dx;mm[1]=tempLine.dy;mm[2]=tempLine.dz;
		mm[3]=1;
		tempLine.dx=(matriZ.A[0][0]*mm[0]+matriZ.A[1][0]*mm[1]+matriZ.A[2][0]*mm[2]+matriZ.A[3][0]*mm[3]);
		tempLine.dy=(matriZ.A[0][1]*mm[0]+matriZ.A[1][1]*mm[1]+matriZ.A[2][1]*mm[2]+matriZ.A[3][1]*mm[3]);
		tempLine.dz=(matriZ.A[0][2]*mm[0]+matriZ.A[1][2]*mm[1]+matriZ.A[2][2]*mm[2]+matriZ.A[3][2]*mm[3]);
		m_vecVECTOR3D.push_back(tempLine);
	}
	m[0]=tempLine.dx;m[1]=tempLine.dy;m[2]=tempLine.dz;
	m[3]=1;
	tempLine.dx=(matriX.A[0][0]*m[0]+matriX.A[1][0]*m[1]+matriX.A[2][0]*m[2]+matriX.A[3][0]*m[3]);
	tempLine.dy=(matriX.A[0][1]*m[0]+matriX.A[1][1]*m[1]+matriX.A[2][1]*m[2]+matriX.A[3][1]*m[3]);
	tempLine.dz=(matriX.A[0][2]*m[0]+matriX.A[1][2]*m[1]+matriX.A[2][2]*m[2]+matriX.A[3][2]*m[3]);
	for (int i=0;i<(360/angB);i++)
	{
		double mm[4];
		mm[0]=tempLine.dx;mm[1]=tempLine.dy;mm[2]=tempLine.dz;
		mm[3]=1;
		tempLine.dx=(matriZ.A[0][0]*mm[0]+matriZ.A[1][0]*mm[1]+matriZ.A[2][0]*mm[2]+matriZ.A[3][0]*mm[3]);
		tempLine.dy=(matriZ.A[0][1]*mm[0]+matriZ.A[1][1]*mm[1]+matriZ.A[2][1]*mm[2]+matriZ.A[3][1]*mm[3]);
		tempLine.dz=(matriZ.A[0][2]*mm[0]+matriZ.A[1][2]*mm[1]+matriZ.A[2][2]*mm[2]+matriZ.A[3][2]*mm[3]);
		m_vecVECTOR3D.push_back(tempLine);
	}
	for (int i=0;i<m_vecVECTOR3D.size();i++)
	{
		double mm[4];
		mm[0]=m_vecVECTOR3D[i].dx;mm[1]=m_vecVECTOR3D[i].dy;mm[2]=m_vecVECTOR3D[i].dz;
		mm[3]=1;
		m_vecVECTOR3D[i].dx=(currentB.A[0][0]*mm[0]+currentB.A[1][0]*mm[1]+currentB.A[2][0]*mm[2]+currentB.A[3][0]*mm[3]);
		m_vecVECTOR3D[i].dy=(currentB.A[0][1]*mm[0]+currentB.A[1][1]*mm[1]+currentB.A[2][1]*mm[2]+currentB.A[3][1]*mm[3]);
		m_vecVECTOR3D[i].dz=(currentB.A[0][2]*mm[0]+currentB.A[1][2]*mm[1]+currentB.A[2][2]*mm[2]+currentB.A[3][2]*mm[3]);
	}
}
void Ccalcubase::create_cone_new(vector<VECTOR3D>& m_vecVECTOR3D,VECTOR3D LineNorZ,VECTOR3D LineNor,POINT3D orig,double angA,int numB)//建立扇形可视锥
{
	m_vecVECTOR3D.clear();
	VECTOR3D tempLine;
	MATRIX3D matriZ,currentB,currentM;
	create_current_Frenet(currentM,LineNorZ,orig);
	create_BACK_Frenet(currentB,LineNorZ,orig);
	create_current_matriZ(matriZ,angA);
	double mm[4];
	mm[0]=LineNor.dx;mm[1]=LineNor.dy;mm[2]=LineNor.dz;
	mm[3]=1;
	LineNor.dx=(currentM.A[0][0]*mm[0]+currentM.A[1][0]*mm[1]+currentM.A[2][0]*mm[2]+currentM.A[3][0]*mm[3]);
	LineNor.dy=(currentM.A[0][1]*mm[0]+currentM.A[1][1]*mm[1]+currentM.A[2][1]*mm[2]+currentM.A[3][1]*mm[3]);
	LineNor.dz=(currentM.A[0][2]*mm[0]+currentM.A[1][2]*mm[1]+currentM.A[2][2]*mm[2]+currentM.A[3][2]*mm[3]);
	mm[0]=LineNor.dx;mm[1]=LineNor.dy;mm[2]=LineNor.dz;
	mm[3]=1;
	tempLine.dx=(matriZ.A[0][0]*mm[0]+matriZ.A[1][0]*mm[1]+matriZ.A[2][0]*mm[2]+matriZ.A[3][0]*mm[3]);
	tempLine.dy=(matriZ.A[0][1]*mm[0]+matriZ.A[1][1]*mm[1]+matriZ.A[2][1]*mm[2]+matriZ.A[3][1]*mm[3]);
	tempLine.dz=(matriZ.A[0][2]*mm[0]+matriZ.A[1][2]*mm[1]+matriZ.A[2][2]*mm[2]+matriZ.A[3][2]*mm[3]);
	tempLine.Normalize();
	create_current_matriZ(matriZ,5);
	for (int i=0;i<numB;i++)
	{
		double mm[4];
		mm[0]=tempLine.dx;mm[1]=tempLine.dy;mm[2]=tempLine.dz;
		mm[3]=1;
		tempLine.dx=(matriZ.A[0][0]*mm[0]+matriZ.A[1][0]*mm[1]+matriZ.A[2][0]*mm[2]+matriZ.A[3][0]*mm[3]);
		tempLine.dy=(matriZ.A[0][1]*mm[0]+matriZ.A[1][1]*mm[1]+matriZ.A[2][1]*mm[2]+matriZ.A[3][1]*mm[3]);
		tempLine.dz=(matriZ.A[0][2]*mm[0]+matriZ.A[1][2]*mm[1]+matriZ.A[2][2]*mm[2]+matriZ.A[3][2]*mm[3]);
		tempLine.Normalize();
		m_vecVECTOR3D.push_back(tempLine);
	}
	for (int i=0;i<m_vecVECTOR3D.size();i++)
	{
		double mm[4];
		mm[0]=m_vecVECTOR3D[i].dx;mm[1]=m_vecVECTOR3D[i].dy;mm[2]=m_vecVECTOR3D[i].dz;
		mm[3]=1;
		m_vecVECTOR3D[i].dx=(currentB.A[0][0]*mm[0]+currentB.A[1][0]*mm[1]+currentB.A[2][0]*mm[2]+currentB.A[3][0]*mm[3]);
		m_vecVECTOR3D[i].dy=(currentB.A[0][1]*mm[0]+currentB.A[1][1]*mm[1]+currentB.A[2][1]*mm[2]+currentB.A[3][1]*mm[3]);
		m_vecVECTOR3D[i].dz=(currentB.A[0][2]*mm[0]+currentB.A[1][2]*mm[1]+currentB.A[2][2]*mm[2]+currentB.A[3][2]*mm[3]);
	}
}
void Ccalcubase::cone_delel(vector<double>& vec_dis)
{
	vector<double>::iterator ite;
	sort(vec_dis.begin(),vec_dis.end());
	//double cent_dis=0;//中位数
	//if (vec_dis.size()%2==0)
	//{
	//	cent_dis=vec_dis[vec_dis.size()*0.5-1]+vec_dis[vec_dis.size()*0.5];
	//	cent_dis=cent_dis*0.5;
	//}
	//else
	//{
	//	cent_dis=vec_dis[(vec_dis.size()-1)*0.5];
	//}
	//double nomal_dis=0;//中位数标准差
	//for (int i=0;i<vec_dis.size();i++)
	//{
	//	nomal_dis=nomal_dis+(vec_dis[i]-cent_dis)*(vec_dis[i]-cent_dis);
	//}
	//nomal_dis=sqrt(nomal_dis);
	//for (int i=0;i<vec_dis.size();i++)
	//{
	//	if ((vec_dis[i]-cent_dis)>nomal_dis&&(vec_dis[i]-cent_dis)<(-nomal_dis))
	//	{
	//		ite=vec_dis.begin()+i;
	//		vec_dis.erase(ite);
	//		i--;
	//	}
	//}
	vector<double> temp;
	for (int i=0;i<vec_dis.size();i++)
	{
		if (vec_dis[i]<2*vec_dis[0])
		{
			temp.push_back(vec_dis[i]);
		}
		else
		{
			break;
		}
	}
	vec_dis.clear();
	vec_dis=temp;
}
double Ccalcubase::calcuPOINTV(PVERT A)
{
	vec_PFACETTRI vecpFacT;
	vec_VECTOR3D vecVTR;
	VECTOR3D N;
	N = Ccalcubase::CalcuVerNormal(A);
	FindOneRFac(A,vecpFacT);
	FindTwoRing(A,vecpFacT);
	for (int i=0;i<vecpFacT.size();i++)
	{
		VECTOR3D P;
		P.dx=Center(vecpFacT[i]->pHEdge).x-A->x;
		P.dy=Center(vecpFacT[i]->pHEdge).y-A->y;
		P.dz=Center(vecpFacT[i]->pHEdge).z-A->z;
		P.Normalize();
		vecVTR.push_back(P);
	}
	double NV=0;
	for (int i=0;i<vecVTR.size();i++)
	{
		if ((N|vecVTR[i])>NV)
		{
			NV=N|vecVTR[i];
		}
	}
	return NV;
}
double Ccalcubase::calcuPOINTV(PFACETTRI A)
{
	vec_PFACETTRI vecpFacP;
	FindPOneRFAC(A,vecpFacP);
	vec_VECTOR3D nVP;
	POINT3D thePOT=Center(A->pHEdge);
	VECTOR3D N=(VECTOR3D)(*A->m_PFacetNorm);
	double NV=0;
	for (int i=0;i<vecpFacP.size();i++)
	{
		POINT3D temP;
		VECTOR3D temV;
		temP=Center(vecpFacP[i]->pHEdge);
		temV=temP-thePOT;
		temV.Normalize();
		nVP.push_back(temV);
	}
	for (int i=0;i<nVP.size();i++)
	{
		NV=(nVP[i]|N)+NV;
	}
	return NV/nVP.size();
}
double Ccalcubase::calcuNoV(PVERT A)
{
	vec_PFACETTRI vecpFacT;
	vec_VECTOR3D vecVTR;
	VECTOR3D N;
	N = Ccalcubase::CalcuVerNormal(A);
	FindOneRFac(A,vecpFacT);
	//FindTwoRing(A,vecpFacT);
	for (int i=0;i<vecpFacT.size();i++)
	{
		VECTOR3D P;
		P.dx=Center(vecpFacT[i]->pHEdge).x-A->x;
		P.dy=Center(vecpFacT[i]->pHEdge).y-A->y;
		P.dz=Center(vecpFacT[i]->pHEdge).z-A->z;
		P.Normalize();
		vecVTR.push_back(P);
	}
	double NV=0;
	for (int i=0;i<vecVTR.size();i++)
	{
		NV=(vecVTR[i]|N)+NV;
	}
	return NV/vecVTR.size();
	/*for (int i=0;i<vecVTR.size();i++)
	{
		if ((N|vecVTR[i])>NV)
		{
			NV=N|vecVTR[i];
		}
	}
	return NV;*/
}
 double Ccalcubase::calcuNoV_new(PVERT A)
 {
	 vec_PFACETTRI vecpFacT;
	vec_VECTOR3D vecVTR;
	VECTOR3D N;
	N = Ccalcubase::CalcuVerNormal(A);
	FindOneRFac(A,vecpFacT);
	//FindTwoRing(A,vecpFacT);
	for (int i=0;i<vecpFacT.size();i++)
	{
		VECTOR3D P;
		P.dx=Center(vecpFacT[i]->pHEdge).x-A->x;
		P.dy=Center(vecpFacT[i]->pHEdge).y-A->y;
		P.dz=Center(vecpFacT[i]->pHEdge).z-A->z;
		P.Normalize();
		vecVTR.push_back(P);
	}
	double NV=0;
	for (int i=0;i<vecVTR.size();i++)
	{
		if ((N|vecVTR[i])>NV)
		{
			NV=N|vecVTR[i];
		}
	}
	return NV;
 }
BOOL Ccalcubase::FindOneRing(PVERT& pVer,vec_PVERT& vecpVer)
{
	PVERT pVerAdj;
	PHEDGE pHEO;
	it_vec_PVERT it_vecPV;
	it_vecPV = vecpVer.begin();
	pHEO = pVer->pHEdgeOut;
	vecpVer.clear();
	do 
	{
		pVerAdj = pHEO->pVertEnd;
		vecpVer.push_back(pVerAdj);             //点发散半边的终点为其一阶领域点

		if (pHEO->pHEdgePre->pHEdgePair != NULL)   //假如没碰到边界
		{ 
			pHEO = pHEO->pHEdgePre->pHEdgePair; //逆时针寻找，点在点表中逆时针排列
		}
		else                                    //假如碰到边界
		{
			pVerAdj = pHEO->pHEdgeNext->pVertEnd;
			vecpVer.push_back(pVerAdj);             //将逆时针方向边界上的点放入一阶领域
			if(pVer->pHEdgeOut->pHEdgePair != NULL) //假如点的发散半边不是最顺时针的半边
			{
				pHEO = pVer->pHEdgeOut->pHEdgePair->pHEdgeNext;//将点下个判断半边转向起始半边的顺时针侧
				do 
				{
					pVerAdj = pHEO->pVertEnd;
					it_vecPV = vecpVer.begin();
					vecpVer.insert(it_vecPV,pVerAdj);    //在容器的开头插入点，使一阶领域点还是逆时针排列
					if (pHEO ->pHEdgePair !=NULL)
					{
						pHEO = pHEO ->pHEdgePair->pHEdgeNext;
					} 
					else
					{
						pHEO = pVer->pHEdgeOut;
						break;
					}

				} while (pHEO != pVer->pHEdgeOut);
			}
			else                       //假如点的发散半边是最顺时针的半边
			{
				pHEO = pVer->pHEdgeOut;
				break;
			}

		}
	} while (pHEO != pVer->pHEdgeOut);          //向外发散的半边不等于起始的半边，继续循环

	if (vecpVer.size() !=0)                     //容器中有点，搜寻成功
	{
		return TRUE;
	}
	return FALSE;
}
BOOL Ccalcubase::FindOnePH(PVERT A,vec_PHEDGE& vecpH)
{
	vecpH.clear();
	PHEDGE temp;
	temp=A->pHEdgeOut;
	vecpH.push_back(temp);
	int k=0;
	while(1)
	{
	    temp=temp->pHEdgePre->pHEdgePair;
		/*if (temp==NULL)
		{
			break;
		}*/
		if (temp==A->pHEdgeOut)
		{
			break;
		}
		/*if (temp==NULL)
		{
			temp=A->pHEdgeOut;
			while (1)
			{
				if (temp->pHEdgePair==NULL)
				{
					k=100;
					break;
				}
			   temp=temp->pHEdgePair->pHEdgeNext;
			   if (temp==NULL)
			   {
				   k=100;
				   break;
			   }
			   if (temp==A->pHEdgeOut)
			   {
				   break;
			   }
			   else
			   {
				   vecpH.push_back(temp);
			   }
			}
		}
		if (k==100)
		{
			break;
		}*/
		if (temp!=NULL)
		{
		  vecpH.push_back(temp);
		}
	}
	return TRUE;
}
BOOL Ccalcubase::FindOneRFac(PVERT pVer,vec_PFACETTRI& vecpFac)
{
	vecpFac.clear();
	PFACETTRI pFacAdj;
	PHEDGE    pHEO;
	it_vec_PFACETTRI it_vecPF;
	it_vecPF = vecpFac.begin();
	pHEO = pVer->pHEdgeOut;

	do 
	{
		pFacAdj = (PFACETTRI)pHEO->pFacetAdj;
		vecpFac.push_back(pFacAdj);             //点发散半边的相邻面为其一阶领域面

		if (pHEO->pHEdgePre->pHEdgePair != NULL)   //假如没碰到边界
		{ 
			pHEO = pHEO->pHEdgePre->pHEdgePair; //逆时针寻找，点在点表中逆时针排列
		}
		else                                    //假如碰到边界
		{
			if(pVer->pHEdgeOut->pHEdgePair != NULL) //假如点的发散半边不是最顺时针的半边
			{
				pHEO = pVer->pHEdgeOut->pHEdgePair->pHEdgeNext;//将面下个判断半边转向起始半边的顺时针侧
				do 
				{
					pFacAdj = (PFACETTRI)pHEO->pFacetAdj;
					it_vecPF = vecpFac.begin();
					vecpFac.insert(it_vecPF,pFacAdj);    //在容器的开头插入面，使一阶领域面还是逆时针排列
					if (pHEO ->pHEdgePair !=NULL)
					{
						pHEO = pHEO ->pHEdgePair->pHEdgeNext;
					} 
					else
					{
						pHEO = pVer->pHEdgeOut;
					}

				} while (pHEO != pVer->pHEdgeOut);
			}
			else
			{
				pHEO = pVer->pHEdgeOut;
			}

		}

	} while (pHEO != pVer->pHEdgeOut);          //向外发散的半边不等于起始的半边，继续循环

	if (vecpFac.size() !=0)                     //容器中有面，搜寻成功
	{
		return TRUE;
	}
	return FALSE;
}
BOOL Ccalcubase::FindPOneRFAC(PFACETTRI pFac,vec_PFACETTRI& vecpFacP)
{
	/*vecpFacP.clear();
	PFACETTRI B;
	if (A->pHEdge->pHEdgePair!=NULL)
	{
		B=(PFACETTRI)(A->pHEdge->pHEdgePair->pFacetAdj);
		vecpFacP.push_back(B);
	}
	if (A->pHEdge->pHEdgeNext->pHEdgePair!=NULL)
	{
		B=(PFACETTRI)(A->pHEdge->pHEdgeNext->pHEdgePair->pFacetAdj);
		vecpFacP.push_back(B);
	}
	if (A->pHEdge->pHEdgePre->pHEdgePair!=NULL)
	{
		B=(PFACETTRI)(A->pHEdge->pHEdgePre->pHEdgePair->pFacetAdj);
		vecpFacP.push_back(B);
	}
	if (!vecpFacP.empty())
	{
		return TRUE;
	}*/
	pFac->becut=1;
	for(int i=0;i<3;i++)
	{
		vec_PFACETTRI vectemp;
		FindOneRFac(pFac->m_PVerts[i],vectemp);
		for (int j=0;j<vectemp.size();j++)
		{
			if (vectemp[j]->becut==0)
			{
				vecpFacP.push_back(vectemp[j]);
				vectemp[j]->becut=1;
			}
		}
	}
	for (int i=0;i<vecpFacP.size();i++)
	{
		vecpFacP[i]->becut=0;
	}
	return TRUE;
}
BOOL Ccalcubase::FindPOneRFAC_NEW(PFACETTRI pFac,vec_PFACETTRI& vecpFacP)
{
	vecpFacP.clear();
	PFACETTRI B;
	if (pFac->pHEdge->pHEdgePair!=NULL)
	{
		B=(PFACETTRI)(pFac->pHEdge->pHEdgePair->pFacetAdj);
		vecpFacP.push_back(B);
	}
	if (pFac->pHEdge->pHEdgeNext->pHEdgePair!=NULL)
	{
		B=(PFACETTRI)(pFac->pHEdge->pHEdgeNext->pHEdgePair->pFacetAdj);
		vecpFacP.push_back(B);
	}
	if (pFac->pHEdge->pHEdgePre->pHEdgePair!=NULL)
	{
		B=(PFACETTRI)(pFac->pHEdge->pHEdgePre->pHEdgePair->pFacetAdj);
		vecpFacP.push_back(B);
	}
	return TRUE;
}
void Ccalcubase::FindPtwoRFAC(PFACETTRI pFac, vec_PFACETTRI& vecpFacP)
{
	vec_PFACETTRI vectempFac;
	FindPOneRFAC(pFac, vectempFac);
	for (int i = 0; i < vectempFac.size();++i)
	{
		vectempFac[i]->becut = true;
		vecpFacP.push_back(vectempFac[i]);
	}
	for (int i = 0; i < vectempFac.size();++i)
	{
		vec_PFACETTRI vectempFac2;
		FindPOneRFAC(vectempFac[i], vectempFac2);
		for (int k = 0; k < vectempFac2.size();++k)
		{
			if (vectempFac2[k]->becut==false)
			{
				vecpFacP.push_back(vectempFac2[k]);
			}
		}
	}
	for (int i = 0; i < vecpFacP.size();++i)
	{
		vecpFacP[i]->becut = 0;
	}
}
double Ccalcubase::AreaTri(PPOINT3D pPoint1, PPOINT3D pPoint2, PPOINT3D pPoint3)
{
	VECTOR3D vEdge1,vEdge2,vArea;
	double Area;
	vEdge1 = *pPoint2 - *pPoint1;
	vEdge2 = *pPoint3 - *pPoint1;
	vArea  = vEdge1*vEdge2;
	Area   = 0.5*vArea.GetLength();

	return Area;
	//return 0;
}
double Ccalcubase::AreaTri(POINT3D Point1,POINT3D Point2,POINT3D Point3)
{
	VECTOR3D vEdge1,vEdge2,vArea;
	double Area;
	vEdge1 = Point2 - Point1;
	vEdge2 = Point3 - Point1;
	vArea  = vEdge1*vEdge2;
	Area   = 0.5*vArea.GetLength();

	return Area;
	//return 0;
}
double Ccalcubase::AreaTri(PFACETTRI pFacTri)
{
	PVERT pVer1,pVer2,pVer3;
	pVer1 = pFacTri->m_PVerts[0];
	pVer2 = pFacTri->m_PVerts[1];
	pVer3 = pFacTri->m_PVerts[2];
	VECTOR3D vEdge1,vEdge2,vArea;
	double Area;
	vEdge1 = *pVer2 - *pVer1;
	vEdge2 = *pVer3 - *pVer1;
	vArea  = vEdge1*vEdge2;
	Area   = 0.5*vArea.GetLength();

	return Area;
}
VECTOR3D Ccalcubase::CalcuVerNormal(PVERT pVer/*,int FacID*/)
{
	VECTOR3D Normal = VECTOR3D(0,0,0);
	VECTOR3D NorFac,V1,V2;
	vec_PVERT     vecpVer;
	vec_PFACETTRI vecpFac;
	double Area,Angle,AA = 0.0;

	//if(pVer->bSharpVer == FALSE)                    //假如点不是锐边点
	//{
	FindOneRing(pVer,vecpVer);
	FindOneRFac(pVer,vecpFac);
	int Nver = vecpVer.size();
	int Nfac = vecpFac.size();
	int nNum = 0;
	if (Nfac == Nver) nNum = Nver;
	else              nNum = Nfac;

	for (int i=0;i<nNum;i++)
	{
		Area  = AreaTri(pVer,vecpVer[(i+1)%Nver],vecpVer[i]);//计算相邻面片的面积
		V1    = *vecpVer[i] - *pVer;
		V2    = *vecpVer[(i+1)%Nver] - *pVer;
		Angle = _AngleBetween3D(V1,V2);                     //计算点发散出去向量的夹角
		NorFac= *(vecpFac[i]->m_PFacetNorm);
		Normal= Normal + NorFac*Area*Angle;
		AA    = AA + Area*Angle;
	}

	Normal = Normal/AA;                                    //法矢加权N = W(NorFac*Area*Angle)/W(Area*Angle)
	Normal = Normal.GetNormal();                           //单位化
	//}
	//else
	//{

	//}
	return Normal;
}
VECTOR3D Ccalcubase::CalcuPntInFacNormal(POINT3D Pt,PFACETTRI Fac)
{
	VECTOR3D Vp,V0,V1,V2;
	POINT3D  P0,P1,P2;
	double   Af,A0,A1,A2;

	V0 = Fac->m_PVerts[0]->Normal;
	V1 = Fac->m_PVerts[1]->Normal;
	V2 = Fac->m_PVerts[2]->Normal;
	P0 = (POINT3D)*(Fac->m_PVerts[0]);
	P1 = (POINT3D)*(Fac->m_PVerts[1]);
	P2 = (POINT3D)*(Fac->m_PVerts[2]);

	Af = AreaTri(Fac);
	A0 = AreaTri(Pt,P1,P2);
	A1 = AreaTri(Pt,P2,P0);
	A2 = AreaTri(Pt,P1,P0);

	Vp = (V0*A0 + V1*A1 + V2*A2)/Af;
	Vp.Normalize();
	return Vp;
}
VECTOR3D Ccalcubase::CalcuPntInLineNormal(POINT3D Pt,PHEDGE theLine)
{
	VECTOR3D theVP,thL1,theL2,theL;
	POINT3D P1,P2;
	double L,L1,L2;
	P1=*theLine->pVertEnd;
	P2=*theLine->pHEdgePair->pVertEnd;
	theL=P1-P2;
	thL1=Pt-P1;
	theL2=Pt-P2;
	L=theL.GetLength();L1=thL1.GetLength();L2=theL2.GetLength();
	theVP.dx=L2*theLine->pVertEnd->Normal.dx +L1*theLine->pHEdgePair->pVertEnd->Normal.dx;
	theVP.dy=L2*theLine->pVertEnd->Normal.dy +L1*theLine->pHEdgePair->pVertEnd->Normal.dy;
	theVP.dz=L2*theLine->pVertEnd->Normal.dz +L1*theLine->pHEdgePair->pVertEnd->Normal.dz;
	theVP.dx=theVP.dx/L;theVP.dy=theVP.dy/L;theVP.dz=theVP.dz/L;
	theVP.Normalize();
	return theVP;
}
//计算顶点一阶领域锐角三角形的混合面积（三角以外心为界分割，与顶点[pPoint1]相邻的部分的面积）
double Ccalcubase::AreaTriMix(PPOINT3D pPoint1, PPOINT3D pPoint2, PPOINT3D pPoint3)   //逆时针排列
{
	double AreaTM;
	VECTOR3D vEdge1,vEdge2,vEdge3;
	vEdge1 = *pPoint2 - *pPoint1;
	vEdge2 = *pPoint3 - *pPoint2;
	vEdge3 = *pPoint1 - *pPoint3;
	double dV1 = vEdge1|vEdge1;                                  //a^2
	double dV3 = vEdge3|vEdge3;
	double dCotA2,dCotA3;
	VECTOR3D vSin2,vSin3;
	vSin2 = (VECTOR3D(0,0,0) -vEdge1)*vEdge2;                         //a*b
	vSin3 = (VECTOR3D(0,0,0) -vEdge2)*vEdge3;

	dCotA2 = ((VECTOR3D(0,0,0) - vEdge1)|vEdge2)/vSin2.GetLength();      //cotA = (a|b)/|a*b|
	dCotA3 = ((VECTOR3D(0,0,0) - vEdge2)|vEdge3)/vSin3.GetLength();

	AreaTM = 0.125*(dV3*dCotA2 + dV1*dCotA3);     //Am = 0.125*(c^2*cotA2 + a^2*cotA3)

	return AreaTM;
}
//混合面积总和
double Ccalcubase::AreaTriMixSum(PVERT pVer)
{
	double AreaTMS=0,AreaM=0;
	VECTOR3D vE1,vE2,vE3;
	double dE1,dE2,dE3,MaxE;               //三角形三条边的长度
	double deta;                             //deta = c^2-(a^2+b^2)
	vec_PVERT     vecpVer;
	vec_PFACETTRI vecpFac;
	//BOOL isBoundV;                         //判断是否为边界点

	BOOL IsOne = FindOneRing(pVer,vecpVer);       //搜寻一阶领域点
	if (IsOne == FALSE) return FALSE;
	FindOneRFac(pVer,vecpFac);

	int nFac = vecpFac.size();
	int nVer = vecpVer.size(); //一阶领域点在容器中逆时针存储
	int nNum = 0;
	if (nFac == nVer) nNum = nVer;//判断是否为边界点
	else              nNum = nFac;
	for (int i=0;i<nNum;i++)
	{
		vE1 = *vecpVer[i] - *pVer;
		vE2 = *vecpVer[(i+1)%nVer] - *vecpVer[i];
		vE3 = *pVer - *vecpVer[(i+1)%nVer];
		dE1 = vE1.GetLength();
		dE2 = vE2.GetLength();
		dE3 = vE3.GetLength();
		MaxE= max(max(dE1,dE2),dE3);

		if (MaxE == dE1)                            //判断三角形是否为锐角三角形
		{	deta = dE1*dE1-dE2*dE2-dE3*dE3;
		} 
		else
		{
			if (MaxE == dE2)
			{	deta = dE2*dE2-dE1*dE1-dE3*dE3;
			} 
			else
			{	deta = dE3*dE3-dE1*dE1-dE2*dE2;
			}
		}

		if (deta<0)
		{   AreaM = AreaTriMix(pVer,vecpVer[i],vecpVer[(i+1)%nVer]);
		} 
		else
		{
			if (MaxE == dE2)
			{	AreaM = 0.5*AreaTri(pVer,vecpVer[i],vecpVer[(i+1)%nVer]);//计算公式参考文献
			} 
			else
			{	AreaM = 0.25*AreaTri(pVer,vecpVer[i],vecpVer[(i+1)%nVer]);
			}
		}

		AreaTMS = AreaTMS + AreaM;
	}
	return AreaTMS;
}
//计算顶点的平均曲率
double Ccalcubase::CalcuVerMeanCurvature(PVERT pVer)
{
	double VerMeanCur=0;
	vec_PVERT     vecpVer;
	vec_PFACETTRI vecpFac;
	BOOL IsOne = FindOneRing(pVer,vecpVer);       //搜寻一阶领域点
	if (IsOne == FALSE) return FALSE;
	FindOneRFac(pVer,vecpFac);
	
	double Am = AreaTriMixSum(pVer);
	VECTOR3D vNor =pVer->Normal;

	VECTOR3D Vij,V1ji,V1jj,V2ji,V2jj;
	double dCotA1,dCotA2;
	VECTOR3D vSin1,vSin2;

	double dVN;

	int nFac = vecpFac.size();
	int nVer = vecpVer.size(); //一阶领域点在容器中逆时针存储
	int nNum = 0;
	if (nFac == nVer) nNum = nVer;//判断是否为边界点
	else              nNum = nFac;
	for (int i=0;i<nNum;i++)
	{
		Vij  = *vecpVer[i] - *pVer;
		V1ji = *pVer       - *vecpVer[(i-1+nVer)%nVer];      //***********当点为边界点的时候该离散方法就无法使用了！！！！！
		V1jj = *vecpVer[i] - *vecpVer[(i-1+nVer)%nVer];
		V2ji = *pVer       - *vecpVer[(i+1)%nVer];
		V2jj = *vecpVer[i] - *vecpVer[(i+1)%nVer];

		vSin1 = V1jj*V1ji;                         //a*b
		vSin2 = V2ji*V2jj;
		dCotA1 = (V1jj|V1ji)/vSin1.GetLength();      //cotA = (a|b)/|a*b|
		dCotA2 = (V2ji|V2jj)/vSin2.GetLength();

		dVN = Vij|vNor;
		VerMeanCur += (dCotA1+dCotA2)*dVN;
	}

	VerMeanCur = 0.25*VerMeanCur/Am;
	return VerMeanCur;
}
//计算顶点的高斯曲率
double Ccalcubase::CalcuVerGausCurvature(PVERT pVer)
{
	double VerGausCur=0;
	vec_PVERT vecpVer;
	BOOL IsOne = FindOneRing(pVer,vecpVer);       //搜寻一阶领域点
	if (IsOne == FALSE) return FALSE;
	double Am = AreaTriMixSum(pVer);

	VECTOR3D Vij1,Vij2;
	double dA,dAM=0;

	int nNum = vecpVer.size();                //一阶领域点在容器中逆时针存储
	for (int i=0;i<nNum;i++)
	{
		Vij1 = *vecpVer[i] - *pVer;
		Vij2 = *vecpVer[(i+1)%nNum] - *pVer;
		dA   = _AngleBetween3D(Vij1,Vij2);
		dAM +=dA;
	}
	VerGausCur = (2*PI - dAM)/Am;

	return VerGausCur;
}
//计算顶点的主曲率
void Ccalcubase::CalcuVerPrinCurvature(PVERT pVer,double& K1,double& K2)
{
	double KH,KG;
	KH = CalcuVerMeanCurvature(pVer);
	KG = CalcuVerGausCurvature(pVer);

	K1 = KH + sqrt(fabs(KH*KH - KG));
	K2 = KH - sqrt(fabs(KH*KH - KG));

}
BOOL Ccalcubase::FindRHEd(PHEDGE A,vec_PFACETTRI& m_vecS)
{
	vec_PFACETTRI vecpFac;
	FindOneRFac(A->pVertEnd,vecpFac);
	for (int i=0;i<vecpFac.size();i++)
	{
		vecpFac[i]->bStatus=1;
	}
	m_vecS=vecpFac;
	vecpFac.clear();
	FindOneRFac(A->pHEdgeNext->pVertEnd,vecpFac);
	for (int i=0;i<vecpFac.size();i++)
	{
		if (vecpFac[i]->bStatus==0)
		{
			vecpFac[i]->bStatus=1;
			m_vecS.push_back(vecpFac[i]);
		}
	}
	vecpFac.clear();
	FindOneRFac(A->pHEdgeNext->pHEdgeNext->pVertEnd,vecpFac);
	for (int i=0;i<vecpFac.size();i++)
	{
		if (vecpFac[i]->bStatus==0)
		{
			vecpFac[i]->bStatus=1;
			m_vecS.push_back(vecpFac[i]);
		}
	}
	vecpFac.clear();
	FindOneRFac(A->pHEdgePair->pHEdgeNext->pVertEnd,vecpFac);
	for (int i=0;i<vecpFac.size();i++)
	{
		if (vecpFac[i]->bStatus==0)
		{
			vecpFac[i]->bStatus=1;
			m_vecS.push_back(vecpFac[i]);
		}
	}
	vecpFac.clear();
	if (m_vecS.size()!=0)
	{
		return TRUE;
	}
	else
	{
		return FALSE;
	}
}
//计算网格点的二阶领域面
BOOL Ccalcubase::FindTwoRing(PVERT pVer,vec_PFACETTRI& vecpFacT)
{
	vec_PVERT  vecpVer;
	PFACETTRI pFacAdj;
	PVERT      pVerAdjOne,pVerAdj;
	set<PFACETTRI> setpFacT;
	pair<set<PFACETTRI>::iterator,bool> pair_insertResult;
	it_vec_PFACETTRI it_vecPV;
	it_vecPV = vecpFacT.begin();

	BOOL IsOne = FindOneRing(pVer,vecpVer);       //搜寻一阶领域点
	if (IsOne == FALSE)
	{
		return FALSE;
	}

	int nNum = vecpVer.size();
	for (int i=0;i<nNum;i++)
	{
		pVerAdjOne = vecpVer[i];
		PHEDGE pHEO;
		pHEO = pVerAdjOne->pHEdgeOut;
		do 
		{	
			pFacAdj = (PFACETTRI)pHEO->pFacetAdj;
			if (pHEO->pHEdgePair == NULL)
			{
				//pVerAdjOne->bStatus= TRUE;
				pFacAdj->bStatus=TRUE;
				break;
			}
			pHEO = pHEO->pHEdgePair->pHEdgeNext;   //顺时针搜索

		} while (pHEO != pVerAdjOne->pHEdgeOut);


		if (pFacAdj->bStatus==TRUE)
		{
			pair_insertResult = setpFacT.insert(pFacAdj);
			if (pair_insertResult.second)
			{
				vecpFacT.push_back(pFacAdj);
			}
		} 
		else
		{
			do 
			{
				pVerAdj = pHEO->pVertEnd;
				pFacAdj = (PFACETTRI)pHEO->pFacetAdj;
				if (pFacAdj != pVer->pHEdgeOut->pFacetAdj)
				{
					pair_insertResult = setpFacT.insert(pFacAdj);
					if (pair_insertResult.second)
					{
						vecpFacT.push_back(pFacAdj);
					}

				}
				pHEO = pHEO->pHEdgePre->pHEdgePair; //逆时针寻找，点在点表中逆时针排列

			} while (pHEO != pVerAdjOne->pHEdgeOut);          //向外发散的半边不等于起始的半边，继续循环
		}
	}

	if (vecpFacT.size() !=0)                     //容器中有点，搜寻成功
	{
		return TRUE;
	}
	return FALSE;
}
//计算网格点的二阶领域点
BOOL Ccalcubase::FindTwoRing(PVERT pVer,vec_PVERT& vecpVerT)
{
	vec_PVERT  vecpVer;
	PVERT      pVerAdjOne,pVerAdj;
	set<PVERT> setpVerT;
	pair<set<PVERT>::iterator,bool> pair_insertResult;
	it_vec_PVERT it_vecPV;
	it_vecPV = vecpVerT.begin();
	vecpVerT.clear();
	BOOL IsOne = FindOneRing(pVer,vecpVer);       //搜寻一阶领域点
	if (IsOne == FALSE)
	{
		return FALSE;
	}

	int nNum = vecpVer.size();
	for (int i=0;i<nNum;i++)
	{
		pVerAdjOne = vecpVer[i];
		PHEDGE pHEO;
		pHEO = pVerAdjOne->pHEdgeOut;
		do 
		{	
			pHEO = pHEO->pHEdgePair->pHEdgeNext;   //顺时针搜索

		} while (pHEO != pVerAdjOne->pHEdgeOut);


		if (pVerAdjOne->bSharpVer == TRUE)
		{
			pair_insertResult = setpVerT.insert(pVerAdjOne);
			if (pair_insertResult.second)
			{
				vecpVerT.push_back(pVerAdjOne);
			}
		} 
		else
		{
			do 
			{
				pVerAdj = pHEO->pVertEnd;
				if (pVerAdj != pVer)
				{
					pair_insertResult = setpVerT.insert(pVerAdj);
					if (pair_insertResult.second)
					{
						vecpVerT.push_back(pVerAdj);
					}

				}
				pHEO = pHEO->pHEdgePre->pHEdgePair; //逆时针寻找，点在点表中逆时针排列

			} while (pHEO != pVerAdjOne->pHEdgeOut);          //向外发散的半边不等于起始的半边，继续循环
		}
	}

	if (vecpVerT.size() !=0)                     //容器中有点，搜寻成功
	{
		return TRUE;
	}
	return FALSE;
}
BOOL Ccalcubase::FindTwoRing_new(PVERT pVer,vec_PVERT& vecpVerT)
{
	vecpVerT.clear();
	vec_PVERT  vecpVer;
	FindOneRing(pVer,vecpVer);
	set<PVERT> thetemp;
	set<PVERT>::iterator tehIT;
	for (int i=0;i<vecpVer.size();i++)
	{
	   vec_PVERT vec_temp;
	   FindOneRing(vecpVer[i],vec_temp);
	   for (int j=0;j<vec_temp.size();j++)
	   {
		   thetemp.insert(vec_temp[j]);
	   }
	}
	for (tehIT=thetemp.begin();tehIT!=thetemp.end();tehIT++)
	{
		PVERT tempV=*tehIT;
		vecpVerT.push_back(tempV);
	}
	return vecpVerT.empty();
}
BOOL Ccalcubase::findThreeRing(PVERT pVer,vec_PVERT& vecpVerR)
{
	vec_PVERT vecpVerT;
	PVERT      pVerAdjOne,pVerAdj;
	FindTwoRing(pVer,vecpVerT);
	set<PVERT> setpVerR;
	pair<set<PVERT>::iterator,bool> pair_insertResult;
	it_vec_PVERT it_vecPV;
	it_vecPV = vecpVerR.begin();
	vecpVerR.clear();
	int nNum = vecpVerT.size();
	for (int i=0;i<nNum;i++)
	{
		pVerAdjOne = vecpVerT[i];
		PHEDGE pHEO;
		pHEO = pVerAdjOne->pHEdgeOut;
		do 
		{	
			if (pHEO->pHEdgePair == NULL)
			{
				pVerAdjOne->bSharpVer = TRUE;
				break;
			}
			pHEO = pHEO->pHEdgePair->pHEdgeNext;   //顺时针搜索

		} while (pHEO != pVerAdjOne->pHEdgeOut);
		////////////
		if (pVerAdjOne->bSharpVer == TRUE)
		{
			pair_insertResult = setpVerR.insert(pVerAdjOne);
			if (pair_insertResult.second)
			{
				vecpVerR.push_back(pVerAdjOne);
			}
		} 
		else
		{
			do 
			{
				pVerAdj = pHEO->pVertEnd;
				if (pVerAdj != pVer)
				{
					pair_insertResult = setpVerR.insert(pVerAdj);
					if (pair_insertResult.second)
					{
						vecpVerR.push_back(pVerAdj);
					}

				}
				pHEO = pHEO->pHEdgePre->pHEdgePair; //逆时针寻找，点在点表中逆时针排列

			} while (pHEO != pVerAdjOne->pHEdgeOut);          //向外发散的半边不等于起始的半边，继续循环
		}
	}

	if (vecpVerR.size() !=0)                     //容器中有点，搜寻成功
	{
		return TRUE;
	}
	return FALSE;
}

POINT3D Ccalcubase::Center(PHEDGE trigle)
{
	POINT3D mindle;
	mindle.x=(trigle->pVertEnd->x+trigle->pHEdgeNext->pVertEnd->x+trigle->pHEdgeNext->pHEdgeNext->pVertEnd->x)/3;
	mindle.y=(trigle->pVertEnd->y+trigle->pHEdgeNext->pVertEnd->y+trigle->pHEdgeNext->pHEdgeNext->pVertEnd->y)/3;
	mindle.z=(trigle->pVertEnd->z+trigle->pHEdgeNext->pVertEnd->z+trigle->pHEdgeNext->pHEdgeNext->pVertEnd->z)/3;
	return mindle;
}
double Ccalcubase::dis(POINT3D A,POINT3D B)
{
	double L=0;
	L=sqrt((A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y)+(A.z-B.z)*(A.z-B.z));
	//L=(A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y)+(A.z-B.z)*(A.z-B.z);
	return L;
}
void Ccalcubase::create_cone_SDF(vector<VECTOR3D>& m_vecVECTOR3D,POINT3D orig,VECTOR3D liner,double angA,double angB)
{
	m_vecVECTOR3D.clear();
	MATRIX3D currentM,matriX,matriZ,currentB;
	VECTOR3D tempLine;
	double m[4];
	create_current_Frenet(currentM,liner,orig);
	create_current_matriX(matriX,angA);
	create_current_matriZ(matriZ,angB);
	create_BACK_Frenet(currentB,liner,orig);
	tempLine.dx=0;tempLine.dy=0;tempLine.dz=-1;
	m[0]=tempLine.dx;m[1]=tempLine.dy;m[2]=tempLine.dz;
	m[3]=1;
	tempLine.dx=(matriX.A[0][0]*m[0]+matriX.A[1][0]*m[1]+matriX.A[2][0]*m[2]+matriX.A[3][0]*m[3]);
	tempLine.dy=(matriX.A[0][1]*m[0]+matriX.A[1][1]*m[1]+matriX.A[2][1]*m[2]+matriX.A[3][1]*m[3]);
	tempLine.dz=(matriX.A[0][2]*m[0]+matriX.A[1][2]*m[1]+matriX.A[2][2]*m[2]+matriX.A[3][2]*m[3]);
	for (int i=0;i<(360/angB);i++)
	{
		double mm[4];
		mm[0]=tempLine.dx;mm[1]=tempLine.dy;mm[2]=tempLine.dz;
		mm[3]=1;
		tempLine.dx=(matriZ.A[0][0]*mm[0]+matriZ.A[1][0]*mm[1]+matriZ.A[2][0]*mm[2]+matriZ.A[3][0]*mm[3]);
		tempLine.dy=(matriZ.A[0][1]*mm[0]+matriZ.A[1][1]*mm[1]+matriZ.A[2][1]*mm[2]+matriZ.A[3][1]*mm[3]);
		tempLine.dz=(matriZ.A[0][2]*mm[0]+matriZ.A[1][2]*mm[1]+matriZ.A[2][2]*mm[2]+matriZ.A[3][2]*mm[3]);
		m_vecVECTOR3D.push_back(tempLine);
	}
	for (int i=0;i<m_vecVECTOR3D.size();i++)
	{
		double mm[4];
		mm[0]=m_vecVECTOR3D[i].dx;mm[1]=m_vecVECTOR3D[i].dy;mm[2]=m_vecVECTOR3D[i].dz;
		mm[3]=1;
		m_vecVECTOR3D[i].dx=(currentB.A[0][0]*mm[0]+currentB.A[1][0]*mm[1]+currentB.A[2][0]*m[2]+currentB.A[3][0]*mm[3]);
		m_vecVECTOR3D[i].dy=(currentB.A[0][1]*mm[0]+currentB.A[1][1]*mm[1]+currentB.A[2][1]*m[2]+currentB.A[3][1]*mm[3]);
		m_vecVECTOR3D[i].dz=(currentB.A[0][2]*mm[0]+currentB.A[1][2]*mm[1]+currentB.A[2][2]*m[2]+currentB.A[3][2]*mm[3]);
	}
}
//////////回溯法遍历

////////////////////////////////////SDF聚类的两个公式
int Ccalcubase::SDF_center(vector<double> &vec_center,double xi)
{
	int CENT=0;
	double tempDIS;
	tempDIS=(vec_center[0]-xi)*(vec_center[0]-xi);
	for (int i=1;i<vec_center.size();i++)
	{
		if (tempDIS>=(vec_center[i]-xi)*(vec_center[i]-xi))
		{
			tempDIS=(vec_center[i]-xi)*(vec_center[i]-xi);
			CENT=i;
		}
	}
	return CENT;
}
double Ccalcubase::SDF_new_center(vector<double> &vec_c)
{
	double distemp1=0;
	double CENT=0;
	for (int j=0;j<vec_c.size();j++)
	{
		distemp1=distemp1+(vec_c[j]-vec_c[0])*(vec_c[j]-vec_c[0]);
	}
	for (int i=1;i<vec_c.size();i++)
	{
		double distemp=0;
		for (int j=0;j<vec_c.size();j++)
		{
			distemp=distemp+(vec_c[j]-vec_c[i])*(vec_c[j]-vec_c[i]);
		}
		if (distemp<distemp1)
		{
			distemp1=distemp;
			CENT=vec_c[i];
		}
	}
	return CENT;
}
void Ccalcubase::SDF_nomal(vector<double>& vec_SDF,map<double,PFACETTRI>& map_SDF)
{
	/*map<double,PFACETTRI> map_temp;
	double disd=0;
	double ai=4;
	ai=log((double)5);
	disd=vec_SDF.back()-vec_SDF[0];
	for (int i=0;i<vec_SDF.size();i++)
	{
		PFACETTRI theface;
		theface=map_SDF[vec_SDF[i]];
		vec_SDF[i]=log((vec_SDF[i]-vec_SDF[0])*4+1)/ai;
		map_temp.insert(pair<double,PFACETTRI>(vec_SDF[i],theface));
	}
	map_SDF.clear();
	map_SDF=map_temp;*/
	map<double,PFACETTRI> map_temp;
	double disd=0;
	double ai=40;
	ai=log((double)10);
	disd=vec_SDF.back()-vec_SDF[0];
	for (int i=0;i<vec_SDF.size();i++)
	{
		PFACETTRI theface;
		theface=map_SDF[vec_SDF[i]];
		vec_SDF[i]=log(vec_SDF[i])/ai;
		map_temp.insert(pair<double,PFACETTRI>(vec_SDF[i],theface));
	}
	map_SDF.clear();
	map_SDF=map_temp;
}
void Ccalcubase::SDF_nomal_old(vector<double>& vec_SDF,map<double,PFACETTRI>& map_SDF)//SDF归一化
{
	map<double,PFACETTRI> map_temp;
	double disd=0;
	double ai=4;
	ai=log((double)5);
	disd=vec_SDF.back()-vec_SDF[0];
	for (int i=0;i<vec_SDF.size();i++)
	{
		PFACETTRI theface;
		theface=map_SDF[vec_SDF[i]];
		vec_SDF[i]=log((vec_SDF[i]-vec_SDF[0])*4+1)/ai;
		map_temp.insert(pair<double,PFACETTRI>(vec_SDF[i],theface));
	}
	map_SDF.clear();
	map_SDF=map_temp;
}
void Ccalcubase:: SDF_nomal(vector<double>& vec_SDF,map<double,PVERT>& map_SDF)
{
	map<double,PVERT> map_temp;
	double disd=0;
	double ai=40;
	ai=log((double)10);
	disd=vec_SDF.back()-vec_SDF[0];
	for (int i=0;i<vec_SDF.size();i++)
	{
		PVERT theface;
		theface=map_SDF[vec_SDF[i]];
		vec_SDF[i]=log(vec_SDF[i])/ai;
		map_temp.insert(pair<double,PVERT>(vec_SDF[i],theface));
	}
	map_SDF.clear();
	map_SDF=map_temp;
}
void Ccalcubase::SDF_delet(vector<double>& vec_SDF)
{
	vector<double>::iterator ite;
	sort(vec_SDF.begin(),vec_SDF.end());
	double cent_SDF=0;//中位数
	if (vec_SDF.size()==30)
	{
		if (vec_SDF.size()%2==0)
		{
			cent_SDF=vec_SDF[vec_SDF.size()*0.5-1]+vec_SDF[vec_SDF.size()*0.5];
			cent_SDF=cent_SDF*0.5;
		}
		else
		{
			cent_SDF=vec_SDF[(vec_SDF.size()-1)*0.5];
		}
	double nomal_SDF=0;//中位数标准差
	for (int i=0;i<vec_SDF.size();i++)
	{
		nomal_SDF=nomal_SDF+(vec_SDF[i]-cent_SDF)*(vec_SDF[i]-cent_SDF);
	}
	nomal_SDF=sqrt(nomal_SDF);
	for (int i=0;i<vec_SDF.size();i++)
	{
		if ((vec_SDF[i]-cent_SDF)>nomal_SDF&&(vec_SDF[i]-cent_SDF)<(-nomal_SDF))
		{
			ite=vec_SDF.begin()+i;
			vec_SDF.erase(ite);
			i--;
		}
	}
	}
}
////////////////////////////////////GMM聚类
void Ccalcubase::GMM_seg(vector<double>& vec_SDF,vector<vector<PFACETTRI>>& Rbox,map<double,PFACETTRI>& map_SDF,int k,vector<double>& vec_u,vector<double>& vec_ct,vector<double>& vec_zmj)
{
	Rbox.clear();
	Rbox.resize(k);
	vector<double> vec_wj;
	for (int j=0;j<k;j++)
	{
		///////初始化E
		double count1,count2;
		count1=count2=0;
		count2=GMM_E(j,vec_u,vec_ct,vec_zmj,vec_SDF,vec_wj);
		while (1)
		{
			GMM_M(vec_SDF,vec_wj,vec_u[j],vec_ct[j],vec_zmj[j]);
			count1=GMM_E(j,vec_u,vec_ct,vec_zmj,vec_SDF,vec_wj);
			if (count1>0.8*count2&&count1<1.2*count2)
			{
				break;
			}
			GMM_M(vec_SDF,vec_wj,vec_u[j],vec_ct[j],vec_zmj[j]);
			count2=GMM_E(j,vec_u,vec_ct,vec_zmj,vec_SDF,vec_wj);
			if (count1>0.8*count2&&count1<1.2*count2)
			{
				break;
			}
		}
	}
	for (int i=0;i<vec_SDF.size();i++ )
	{
		int numK=0;
		numK=GMM_FacetoZi(vec_u,vec_ct,vec_zmj,vec_SDF[i]);
		Rbox[numK].push_back(map_SDF[vec_SDF[i]]);
	}
}
double Ccalcubase::GMM_norm(double x,double u,double ct)
{
	return (1/ sqrt(2 * 3.141593 * ct))* pow(2.718282, -(x - u) * (x - u) / (2 * ct));
}
double Ccalcubase::GMM_E(int zi,vector<double>& vec_u,vector<double>& vec_ct,vector<double>& vec_zmj,vector<double>& vec_SDF,vector<double>& vec_wj)
{
	vec_wj.clear();
	for (int i=0;i<vec_SDF.size();i++)
	{
		double bys=0;
		for (int j=0;j<vec_zmj.size();j++)
		{
			bys=GMM_norm(vec_SDF[i],vec_u[j],vec_ct[j])*vec_zmj[j]+bys;
		}
		double tempWJ=0;
		tempWJ=GMM_norm(vec_SDF[i],vec_u[zi],vec_ct[zi])*vec_zmj[zi]/bys;
		vec_wj.push_back(tempWJ);
	}
	double ai=0;
	for (int i=0;i<vec_wj.size();i++)
	{
		ai=ai+vec_wj[i];
	}
	return ai;
}
void Ccalcubase::GMM_M(vector<double>& vec_SDF,vector<double>& vec_wj,double& u,double& ct,double& mj)
{
	mj=0;
	u=0;
    mj=0;
	double ai=0;
	for (int i=0;i<vec_wj.size();i++)
	{
		ai=ai+vec_wj[i];
		u=u+vec_SDF[i]*vec_wj[i];
	}
	mj=ai/vec_wj.size();
	u=u/ai;
	for (int i=0;i<vec_SDF.size();i++)
	{
		ct=ct+vec_wj[i]*(vec_SDF[i]-u)*(vec_SDF[i]-u);
	}
	ct=ct/ai;
	ct=abs(ct);
}
int Ccalcubase::GMM_FacetoZi(vector<double>& vec_u,vector<double>& vec_ct,vector<double>& vec_zmj,double& theface)
{
	int num=0;
	double temp=0;
	for (int i=0;i<vec_u.size();i++)
	{
		double temp1=0;
		temp1=GMM_norm(theface,vec_u[i],vec_ct[i]);
		if (temp1>temp)
		{
			temp=temp1;
			num=i;
		}
	}
	return num;
}
///////////////////////////////////压缩点云模型分块
//void Ccalcubase::segment_cloud(vec_PVERT& sdf_point,vector<vector<PFACETTRI>>& Rbox,vector<PFACETTRI>& m_vecPFacetTri)
//{
//	KDtree* KD;
//	vector<point> vec_p;
//	vector<const float *> knn;
//	for (int i=0;i<sdf_point.size();i++)
//	{
//		point vi(sdf_point[i]->x,sdf_point[i]->y,sdf_point[i]->z);
//		vec_p.push_back(vi);
//	}
//	float* v0=&vec_p[0][0];
//	KD = new KDtree(v0,sdf_point.size());
//	//////////////////////
//	double tehdis=dis((POINT3D)*sdf_point[0],(POINT3D)*sdf_point[1]);
//	for (int i=2;i<sdf_point.size();i++ )
//	{
//		double tempdis=dis((POINT3D)*sdf_point[i],(POINT3D)*sdf_point[i-1]);
//		if (tempdis<tehdis)
//		{
//			tehdis=tempdis;
//		}
//	}
//	while (1)
//	{
//		queue<point> Q;
//		vector<PFACETTRI> tempface;
//		int k=0;
//		for (int i=0;i<sdf_point.size();i++)
//		{
//			if (sdf_point[i]->bStatus==0)
//			{
//				k++;
//				sdf_point[i]->bStatus=1;
//				Q.push(vec_p[i]);
//				break;
//			}
//		}
//		if (k==0)
//		{
//			break;
//		}
//		while(!Q.empty())
//		{
//			KD->find_pt_in_Range(knn,Q.front(),50000*tehdis);
//			for (int i=0;i<knn.size();i++)
//			{
//				int my_ID=0;
//				my_ID=(knn[i]-&vec_p[0][0])/3;
//				if (sdf_point[my_ID]->bStatus==0)
//				{
//					sdf_point[my_ID]->bStatus=1;
//					point vi(vec_p[my_ID]);
//					tempface.push_back(m_vecPFacetTri[my_ID]);
//					Q.push(vi);
//				}
//			}
//			Q.pop();
//		}
//		Rbox.push_back(tempface);
//	}
//
//}
int Ccalcubase::Center_kmeans_cloud(vector<int>& cloud,vec_PVERT& sdf_point)
{
	double distemp1=0;
	int CENT;
	for (int j=0;j<cloud.size();j++)
	{
		distemp1=distemp1+dis((POINT3D)*sdf_point[cloud[j]],(POINT3D)*sdf_point[cloud[0]]);
	}
	for (int i=1;i<cloud.size();i++)
	{
		double distemp=0;
		for (int j=0;j<cloud.size();j++)
		{
			distemp=distemp+dis((POINT3D)*sdf_point[cloud[i]],(POINT3D)*sdf_point[cloud[i-1]]);
		}
		if (distemp<distemp1)
		{
			CENT=cloud[i];
		}
	}
	return CENT;
}
int Ccalcubase::cloud_center(vector<int>& vec_center,PVERT& xin,vec_PVERT& sdf_point)
{
	int CENT=0;
	double tempDIS;
	tempDIS=dis((POINT3D)*sdf_point[vec_center[0]],(POINT3D)*xin);
	for (int i=1;i<vec_center.size();i++)
	{
		double thedis=dis((POINT3D)*sdf_point[vec_center[i]],(POINT3D)*xin);
		if (tempDIS>=thedis)
		{
			tempDIS=thedis;
			CENT=i;
		}
	}
	return CENT;
}
void Ccalcubase::segment_cloud(vec_PVERT& sdf_point,vector<vector<PFACETTRI>>& Rbox,vector<PFACETTRI>& m_vecPFacetTri)
{
	vector<int> cloud_seg[8];
	vector<int> vec_center;
	vec_center.push_back(1);
	vec_center.push_back(10);
	vec_center.push_back(100);
	vec_center.push_back(200);
	vec_center.push_back(300);
	vec_center.push_back(400);
	vec_center.push_back(500);
	vec_center.push_back(600);

	while(1)
	{
		double num=0;
		//////////先分类
		for (int i=0;i<sdf_point.size();i++)
		{
			int c=cloud_center(vec_center,sdf_point[i],sdf_point);
			cloud_seg[c].push_back(i);
		}
		/////////每个类里找新聚点
		for (int i=0;i<8;i++)
		{
			int ki;
			ki=Center_kmeans_cloud(cloud_seg[i],sdf_point);
			if (dis((POINT3D)*sdf_point[ki],(POINT3D)*sdf_point[vec_center[i]])<1)
			{
				vec_center[i]=ki;
				num++;
			}
			else
			{
				vec_center[i]=ki;
			}
		}
		if (num==8)
		{
			break;
		}
	}
	for (int i=0;i<8;i++)
	{
		vec_PFACETTRI tempface;
		for (int j=0;j<cloud_seg[i].size();j++)
		{
			tempface.push_back(m_vecPFacetTri[cloud_seg[i][j]]);
		}
		Rbox.push_back(tempface);
	}
}
double Ccalcubase::CalPointOnLineCurv(VECTOR3D &Tang, PHEDGE &pHEdge, POINT3D &p)
{
	PVERT  A = pHEdge->pVertEnd;  PVERT B = pHEdge->pHEdgePre->pVertEnd;
	VECTOR3D t1,t2;
	double k1(0),k2(0);
	CalculateCurvByTNT(A,t1,t2,k1,k2);
	VECTOR3D ta = t2;
	double  angleA = _AngleBetween3D(Tang, ta);   // 切方向与最小主方向夹角
	double  ka = k2 * cos(angleA) * cos(angleA) + k1 * sin(angleA) * sin(angleA);

	CalculateCurvByTNT(B,t1,t2,k1,k2);
	VECTOR3D  tb = t2;
	double  angleB = _AngleBetween3D(Tang, tb);
	double  kb = k2 * cos(angleB) * cos(angleB) + k1 * sin(angleB) * sin(angleB);

	double Lab = _DistOf(*A, *B);
	double Lpa = _DistOf(p, *A);
	double t = Lpa / Lab;
	ASSERT(t <= 1.0);
	double k = ka * (1.0 - t) + kb * t;
	return k;
}
double Ccalcubase::CalculatePoint_Curv(PFACETTRI pFacet, POINT3D p, VECTOR3D vec)
{
	double  k[3];
	double  Area[3];
	double a;
	double  AreaFac = FacetTriArea(pFacet, a);
	for (int i = 0; i < 3; i ++)
	{
		VECTOR3D t1,t2;
		double k1,k2;
		CalculateCurvByTNT(pFacet->m_PVerts[i],t1,t2,k1,k2);
		double angle = _AngleBetween3D(vec, t2);
		k[i] = k2 * cos(angle) * cos(angle) + k1 * sin(angle) * sin(angle);
		Area[i] = TriArea(p, *pFacet->m_PVerts[(i+1)%3], *pFacet->m_PVerts[(i+2)%3]);
	}
	double u = Area[0]/AreaFac; double v = Area[1]/AreaFac; double w = Area[2]/AreaFac;
	//ASSERT(IS_ZERO(u+v+w-1.0));
	//double test = u+v+w;
	double pk = k[0] * u + k[1] * v + k[2] * w;
	return pk;
}
void  Ccalcubase::CalculateCurvByTNT(PVERT pVex,VECTOR3D &t1,VECTOR3D &t2,double &K1,double &K2)
{
	// 创建平移矩阵T
	VECTOR3D  VecT = POINT3D(0, 0, 0) - *pVex;  // 平移向量
	MATRIX3D  T = MATRIX3D::CreateTransfMatrix(VecT);
	// 创建旋转矩阵R
	VECTOR3D VexNormal;
	//CalculateVexNor(pstlmodel, pVex, VexNormal);
	//CalculateVexNor(pVex, VexNormal);
	VexNormal=CalcuVerNormal(pVex);
	VECTOR3D  u =  (VexNormal * VECTOR3D(0, 0, 1)).GetNormal();   // 单位轴向量u
	VECTOR3D  v = (VexNormal * u).GetNormal();                 // 两正交单位向量叉积，得v，v为单位向量
	double angle = _AngleBetween3D(VECTOR3D(0, 0, 1), VexNormal);
	MATRIX3D R;             // R初始为单位矩阵
	if(!u.IsZeroLength())    // 若u不为零，说明v点法矢不与z轴平行
	{
		R = MATRIX3D::CreateRotateMatrix(angle, u);
	}
	// 计算v点的2阶邻域点在局部坐标下的坐标值
	vec_PVERT  vecPV;
	//VexTworingVex(pstlmodel, pVex, vecPV);      
	VexTwoNeighVex(pVex, vecPV);
	//VexThreeNeighVex(pVex, vecPV);
	Array2D<double>  equ(6, 6, 0.0);    // 矩阵中元素值初始为0
	Array1D<double>  B(6, 0.0);
	for (unsigned int i = 0; i < vecPV.size(); i ++)
	{
		PVERT pv = vecPV.at(i);
		POINT3D pj = (*pv)*T*R;
		double u = pj.x; double v = pj.y; double h = pj.z;             // 邻点局部坐标值 u,v,h
		// 利用邻点的局部坐标通过最小二乘法拟合二次曲面
		// 解线性方程组
		equ[0][0] += u*u*u*u; equ[0][1] += u*u*u*v; equ[0][2] += u*u*v*v; 
		equ[0][3] += u*u*u;     equ[0][4] += u*u*v;     equ[0][5] += u*u;

		equ[1][0] += u*u*u*v; equ[1][1] += u*u*v*v; equ[1][2] += u*v*v*v; 
		equ[1][3] += u*u*v;     equ[1][4] += u*v*v;     equ[1][5] += u*v;

		equ[2][0] += u*u*v*v; equ[2][1] += u*v*v*v; equ[2][2] += v*v*v*v; 
		equ[2][3] += u*v*v;     equ[2][4] += v*v*v;     equ[2][5] += v*v;

		equ[3][0] += u*u*u;     equ[3][1] += u*u*v;     equ[3][2] += u*v*v;   
		equ[3][3] += u*u;         equ[3][4] += u*v;         equ[3][5] += u;

		equ[4][0] += u*u*v;     equ[4][1] += u*v*v;     equ[4][2] += v*v*v;   
		equ[4][3] += u*v;         equ[4][4] += v*v;         equ[4][5] += v;

		equ[5][0] += u*u;         equ[5][1] += u*v;         equ[5][2] += v*v;     
		equ[5][3] += u;             equ[5][4] += v;             equ[5][5] += 1;

		B[0] += u * u * h;  B[1] += u * v * h;    B[2] += v * v * h; 
		B[3] += u * h;        B[4] += v*h;            B[5] += h;
	}

	Array1D<double>  solu(6);
	LU<double>  LUsolve(equ);
	solu = LUsolve.solve(B);
	double a, b, c, p, q, m;
	a = solu[0]; b = solu[1]; c = solu[2]; p = solu[3]; q = solu[4]; m = solu[5];
	VECTOR3D T1, T2;
	double k1(0),k2(0);
	CalculateCurvature(a, b, c, p, q, T1, T2,k1,k2);
	MATRIX3D inverseR = R.TranMatrix();    // R的转秩矩阵
	t1 = T1 * inverseR;
	t2 = T2 * inverseR;
	K1=k1;
	K2=k2;
}
void  Ccalcubase::VexTwoNeighVex(PVERT pVex, vec_PVERT& vecpVex)
{
	vec_PVERT  vecPV;
	vec_PVERT  vecPV0;
	VexNeighVex(pVex, vecPV);
	for (unsigned int i = 0; i < vecPV.size(); i ++)
	{
		VexNeighVex(vecPV.at(i), vecPV0);
		it_vec_PVERT  iterPV;
		for (unsigned int j = 0; j < vecPV0.size(); j ++)
		{
			iterPV = find(vecpVex.begin(), vecpVex.end(), vecPV0.at(j));
			if (iterPV == vecpVex.end() && vecPV0.at(j) != pVex)      // 顶点不重复，且当前点不等于中心点
			{
				vecpVex.push_back(vecPV0.at(j));
			}	
		}
		vecPV0.clear();
	}
}
void  Ccalcubase::VexNeighVex(PVERT pVex, vec_PVERT& vecpVex)
{
	//PVERT  pv = pVex->pHEdgeOut->pVertEnd;
	PHEDGE pHEdge = pVex->pHEdgeOut;
	do 
	{
		vecpVex.push_back(pHEdge->pVertEnd);
		pHEdge = pHEdge->pHEdgePair->pHEdgeNext;
	} while (pHEdge != pVex->pHEdgeOut);
}
void  Ccalcubase::CalculateCurvature(double a, double b, double c, double p, double q, VECTOR3D& T1, VECTOR3D& T2,double &K1,double &K2)
{
	double ContKG = (4 * a * c - b * b) / (( p * p + q * q + 1) * ( p * p + q * q + 1));
	double ContKH = (c + c * p * p + a + a * q * q - b * p * q) / pow(( p * p + q * q + 1), ( 3.0 / 2.0 ));
	K1 = ContKH + sqrt(ContKH * ContKH - ContKG);
	K2 = ContKH - sqrt(ContKH * ContKH - ContKG);
	//第一基本齐式
	double E = 1 + p * p;
	double F = p * q;
	double G = 1 + q * q;
	//第二基本齐式
	double L, M, N;
	L = 2 * a / sqrt(p * p + q * q + 1);
	M = b / sqrt(p * p + q * q + 1);
	N = 2 * c / sqrt(p * p + q * q + 1);
	//K1 = (L * G - 2 * M * F + N * E) / (2 * (E * G - F * F)) + sqrt(pow( (L * G - 2 * M * F + N * E ), 2) - 
	//	4 * (E * G - F * F) * (L * N - M * M)) / (2 * (E * G - F * F));
	//K2 = (L * G - 2 * M * F + N * E) / (2 * (E * G - F * F)) - sqrt(pow( (L * G - 2 * M * F + N * E), 2) - 
	//	4 * (E * G - F * F) * (L * N - M * M)) / (2 * (E * G - F * F));
	T1 = VECTOR3D(G * M - F * N, (E * N - G * ContKH) / 2 + (E * G - F * F) * sqrt(ContKH * ContKH - ContKG), 0);
	T2 = VECTOR3D((E * N - G * ContKH) / 2 + (E * G - F * F) * sqrt(ContKH * ContKH - ContKG), F * L - E * M, 0);
}
double Ccalcubase::FacetTriArea(PFACETTRI pFacet, double& InnerAngle)
{
	double len[3];
	VECTOR3D vec[3];
	for (int i = 0; i < 3; i++)
	{
		len[i] = _DistOf(*pFacet->m_PVerts[i], *pFacet->m_PVerts[(i + 1) % 3]);
		vec[i] = *pFacet->m_PVerts[i] - *pFacet->m_PVerts[(i + 1) % 3];
	}  
	InnerAngle = _AngleBetween3D(vec[0], vec[2]);
	double FacetArea = 0.5 * len[0] * len[2] * sin(fabs(InnerAngle));
	return FacetArea;
}
double Ccalcubase::TriArea(POINT3D p1, POINT3D p2, POINT3D p3)
{
	double len[3];
	VECTOR3D vec[3];
	vec[0] = p2 - p1;
	vec[2] = p3 - p1;
	len[0] = _DistOf(p1, p2);
	len[2] = _DistOf(p1, p3);
	double InnerAngle = _AngleBetween3D(vec[0], vec[2]);
	double FacetArea = 0.5 * len[0] * len[2] * sin(fabs(InnerAngle));
	if (IS_ZERO(FacetArea))
	{
		return 0.0;
	}
	return FacetArea;
}
void Ccalcubase::CalculatePointNor(PFACETTRI pFacet, POINT3D p, VECTOR3D& pNormal)
{
	VECTOR3D  VexNormal[3];
	double  Area[3];   // 与顶点连线形成的三部分面积
	double a = 0.0;
	double  AreaFac = FacetTriArea(pFacet, a);
	for (int i = 0; i < 3; i++)
	{
		CalculateVexNor(pFacet->m_PVerts[i], VexNormal[i]);
		Area[i] = TriArea(p, *pFacet->m_PVerts[(i + 1) % 3], *pFacet->m_PVerts[(i + 2) % 3]);
	}
	double u = Area[0] / AreaFac; double v = Area[1] / AreaFac; double w = Area[2] / AreaFac;
	//ASSERT(IS_ZERO(u+v+w-1.0));
	//double test = u+v+w;
	pNormal = VexNormal[0] * u + VexNormal[1] * v + VexNormal[2] * w;
	pNormal.Normalize();
}
void  Ccalcubase::Facet_nPNeighFacet(PFACETTRI pFacet, vec_PFACETTRI& vecpFacet)
{
	vec_PVERT  vecPV;
	vec_PFACETTRI  vecPF;
	FacetNeighVex(pFacet, vecPV);
	vecpFacet.push_back(pFacet);
	for (unsigned int i = 0; i < vecPV.size(); i++)
	{
		VexNeighFacet(vecPV.at(i), vecPF);
		for (unsigned int j = 0; j < vecPF.size(); j++)
		{
			it_vec_PFACETTRI iterPF;
			iterPF = find(vecpFacet.begin(), vecpFacet.end(), vecPF.at(j));
			if (iterPF == vecpFacet.end())
			{
				vecpFacet.push_back(vecPF.at(j));
			}
		}
		vecPF.clear();
	}
}
bool Ccalcubase::PointInTri(PFACETTRI pFacet, POINT3D P)
{
	POINT3D A = *pFacet->m_PVerts[0];
	POINT3D B = *pFacet->m_PVerts[1];
	POINT3D C = *pFacet->m_PVerts[2];
	double d1 = SameSide(A, B, C, P);
	double d2 = SameSide(B, C, A, P);
	double d3 = SameSide(C, A, B, P);
	if (IS_ZERO(d1* d2 * d3) || (d1* d2 * d3 < -CAD_ZERO))   // 边上
	{
		return false;
	}
	else if (d1 > CAD_ZERO && d2 > CAD_ZERO && d3 > CAD_ZERO)
	{
		return true;
	}
	//return  SameSide(A, B, C, P) && SameSide(B, C, A, P) && SameSide(C, A, B, P);
}
void Ccalcubase::CalculateNor(PFACETTRI pFacet, POINT3D p, VECTOR3D& pNormal)
{
	PHEDGE pHEdge;
	//double k;
	if (PointInTri(pFacet, p))
	{
		//pNormal = CalPointInTriNorByCore(pstlmodel, pFacet, p);
		pNormal = CalPointInTriNorByCore(pFacet, p);
		//k = CalPointInFacCurv(pstlmodel, pFacet, );
	}
	else if (PointOnLine(pFacet, p, pHEdge))
	{
		CalOnLineNorByFacNor(pHEdge, pNormal);
	}
	else
	{
		PVERT pv = NULL;
		VECTOR3D vec[3];
		for (int i = 0; i < 3; i++)
		{
			vec[i] = *pFacet->m_PVerts[i] - p;
			if (IS_ZERO(vec[i].GetLength()))
			{
				pv = pFacet->m_PVerts[i];
				break;
			}
		}
		//CalculateVexNor(pstlmodel, pv, pNormal);
		CalculateVexNor(pv, pNormal);
	}
	pNormal.Normalize();
}

bool Ccalcubase::IntersectLinePlaneYT(PFACETTRI pFacTri, POINT3D LinePt, VECTOR3D LineNor, POINT3D& InterP)
{
	double t;             //直线参数方程x = LinePt.x+ LineNor.x * t; y = LinePt.y+ LineNor.y * t; z = LinePt.z+ LineNor.z * t
	POINT3D Pt;            //点法式平面
	VECTOR3D FacNor;
	Pt = *(pFacTri->m_PVerts[0]);
	FacNor = *(pFacTri->m_PFacetNorm);

	POINT3D InterP_temp;  //平面交点
	VECTOR3D V01, V12, V20, V0i, V1i, V2i, c0, c1, c2;
	bool  IsInFac = false;    //判断是否在三角形内

	if (!IS_ZERO(LineNor | FacNor))
	{
		double d = (Pt - LinePt) | FacNor;
		double angle = _AngleBetween3D(FacNor, LineNor);
		t = d / cos(angle);
		//t = - ((Pt - LinePt) | FacNor) / (LineNor | FacNor);   //t = ((Pt.x – LinePt.x)*FacNor.x+(Pt.y – LinePt.y)*FacNor.y+(Pt.z – LinePt.z)*FacNor.z) 
		//   /(FacNor.x* LineNor.x+ FacNor.y* LineNor.y+ FacNor.z* LineNor.z)    
		InterP_temp = LinePt + LineNor * t;       //将求得参数t带入直线参数方程求交点

		V01 = *(pFacTri->m_PVerts[1]) - *(pFacTri->m_PVerts[0]);
		V12 = *(pFacTri->m_PVerts[2]) - *(pFacTri->m_PVerts[1]);
		V20 = *(pFacTri->m_PVerts[0]) - *(pFacTri->m_PVerts[2]);
		V0i = InterP_temp - *(pFacTri->m_PVerts[0]);
		V1i = InterP_temp - *(pFacTri->m_PVerts[1]);
		V2i = InterP_temp - *(pFacTri->m_PVerts[2]);

		c0 = V01 * V0i;  c1 = V12 * V1i;  c2 = V20 * V2i;
		if (((c0 | c1) > 0.0) && ((c1 | c2) > 0.0) && ((c0 | c2) > 0.0))
		{
			IsInFac = true;
		}
		if ((IS_ZERO(c0.GetLength()) && ((c1 | c2) > 0.0)) ||
			(IS_ZERO(c1.GetLength()) && ((c0 | c2) > 0.0)) ||
			(IS_ZERO(c2.GetLength()) && ((c0 | c1) > 0.0)))
		{
			IsInFac = true;
		}

		//if (((V01*V0i).dz > 0.0 || IS_ZERO((V01*V0i).dz)) &&
		//	 ((V12*V1i).dz > 0.0 || IS_ZERO((V12*V1i).dz)) &&
		//	 ((V20*V2i).dz > 0.0 || IS_ZERO((V20*V2i).dz))) 
		//	 IsInFac = true;        //若叉乘的符号相同则点在面片上
		//if (((V01*V0i).dz < 0.0 || IS_ZERO((V01*V0i).dz)) &&
		//	 ((V12*V1i).dz < 0.0 || IS_ZERO((V12*V1i).dz)) &&
		//	 ((V20*V2i).dz < 0.0 || IS_ZERO((V20*V2i).dz))) 
		//	 IsInFac = true;

		if (IsInFac)
		{
			InterP = InterP_temp;  //平面上交点为面片上交点
			return true;
		}
	}
	return false;         //无交点
}
bool Ccalcubase::Line_TriIntersect(PFACETTRI pFacet, POINT3D LineP, VECTOR3D LineDir, POINT3D& IntersectP)
{
	VECTOR3D vec = *pFacet->m_PVerts[0] - LineP;
	VECTOR3D vec1 = LineP - POINT3D(0, 0, 0);
	double triD = -(*pFacet->m_PFacetNorm | vec);
	double tempD = *pFacet->m_PFacetNorm | LineDir;
	double tempU = ((*pFacet->m_PFacetNorm) | vec);
	if (IS_ZERO(tempD))   // 直线与面片平行
	{
		return false;
	}
	double t = -tempU / tempD;
	IntersectP = LineP + LineDir * t;
	PHEDGE pHEdge;
	if (PointInTri(pFacet, IntersectP))
	{
		return true;
	}
	else if (PointOnLine(pFacet, IntersectP, pHEdge))
	{
		return true;
	}
	return false;
}
CPlane CPlane::BuildPlane(POINT3D a, POINT3D b, POINT3D c)
{
	CPlane  P;
	P.m_PointIn = a;
	P.m_Normal = ((b - a) * (c - a));
	P.m_Normal.Normalize();
	P.d = P.m_Normal | (a - POINT3D(0, 0, 0));
	return P;
}

CPlane CPlane::BuildPlane(VECTOR3D n, POINT3D p)   // p为平面上一点
{
	CPlane P;
	P.m_Normal = n.GetNormal();
	P.d = n | (p - POINT3D(0, 0, 0));
	P.m_PointIn = p;
	return P;
}

CPlane CPlane::BuildPlane(POINT3D p1, POINT3D p2, VECTOR3D NorF)
{
	CPlane Plane;
	Plane.m_PointIn = p1;
	VECTOR3D vec = p2 - p1;
	Plane.m_Normal = vec * NorF;
	Plane.m_Normal.Normalize();
	Plane.d = Plane.m_Normal | (p1 - POINT3D(0, 0, 0));
	return Plane;
}
bool Ccalcubase::PointOnLine(PFACETTRI pFacet, POINT3D P, PHEDGE& pHEdge)
{
	VECTOR3D v0 = *pFacet->m_PVerts[0] - P;
	VECTOR3D v1 = *pFacet->m_PVerts[1] - P;
	VECTOR3D v2 = *pFacet->m_PVerts[2] - P;
	//ASSERT(v0.GetLength() > 0.0 && v1.GetLength() > 0.0 && v2.GetLength() > 0.0);
	if (IS_ZERO(v0.GetLength()) || IS_ZERO(v1.GetLength()) || IS_ZERO(v2.GetLength()))
	{
		return false;
	}

	VECTOR3D v3 = v0 * v1;
	VECTOR3D v4 = v0 * v2;
	VECTOR3D v5 = v1 * v2;
	if (IS_ZERO(v3.GetLength()))
	{
		pHEdge = pFacet->pHEdge->pHEdgeNext;
		return true;
	}
	else if (IS_ZERO(v4.GetLength()))
	{
		pHEdge = pFacet->pHEdge;
		return true;
	}
	else if (IS_ZERO(v5.GetLength()))
	{
		pHEdge = pFacet->pHEdge->pHEdgePre;
		return true;
	}
	return false;
}
void  Ccalcubase::CalculateVexNor(PVERT  pVex, VECTOR3D& VexNormal)
{
	//VECTOR3D  VexNorN;  // 分子
	//double Dernenner;  // 分母
	vec_PFACETTRI vec_pFacet;
	//VexOneringFacet(pstlmodel, pVex, vec_pFacet);
	VexNeighFacet(pVex, vec_pFacet);
	for (unsigned int i = 0; i < vec_pFacet.size(); i++)   //遍历点的邻接面片
	{
		double a;
		double s = FacetTriArea(vec_pFacet.at(i), a);
		vec_pFacet.at(i)->m_PFacetNorm->Normalize();
		VECTOR3D  temp1 = *(vec_pFacet.at(i)->m_PFacetNorm) * a * s;
		VexNormal += temp1;
	}
	VexNormal.Normalize();
}
void  Ccalcubase::FacetNeighVex(PFACETTRI pFacet, vec_PVERT& vecpVex)
{
	for (int i = 0; i < 3; i++)
	{
		vec_PVERT  temp;
		VexNeighVex(pFacet->m_PVerts[i], temp);
		for (unsigned int j = 0; j < temp.size(); j++)
		{
			it_vec_PVERT  iterPV;
			iterPV = find(vecpVex.begin(), vecpVex.end(), temp.at(j));
			if (iterPV == vecpVex.end() && temp.at(j) != pFacet->m_PVerts[(i + 1) % 3] && temp.at(j) != pFacet->m_PVerts[(i + 2) % 3])
			{
				vecpVex.push_back(temp.at(j));
			}
		}
		temp.clear();
	}
}
double Ccalcubase::SameSide(POINT3D A, POINT3D B, POINT3D C, POINT3D P)
{	
	VECTOR3D v0 = B - A;  VECTOR3D v1 = C - A;  VECTOR3D v2 = P - A;
	v0.Normalize(); v1.Normalize(); v2.Normalize();
	VECTOR3D v3 = v0 * v1;
	VECTOR3D v4 = v0 * v2;
	double d = v3 | v4;
	//return d >=0;
	return d;
}
void  Ccalcubase::VexNeighFacet(PVERT pVex, vec_PFACETTRI& vecpFacet)
{
	PHEDGE pHEdge = pVex->pHEdgeOut;
	do
	{
		vecpFacet.push_back((PFACETTRI)pHEdge->pFacetAdj);
		pHEdge = pHEdge->pHEdgePair->pHEdgeNext;
	} while (pHEdge != pVex->pHEdgeOut);
}
VECTOR3D Ccalcubase::CalPointInTriNorByCore(PFACETTRI pFacet, POINT3D p)
{
	VECTOR3D vec[3];
	for (int i = 0; i < 3; i++)
	{
		//CalculateVexNor(pstlmodel, pFacet->m_PVerts[i],  vec[i]);
		CalculateVexNor(pFacet->m_PVerts[i], vec[i]);
	}
	double a;
	double S = FacetTriArea(pFacet, a);
	double S1 = TriArea(p, *pFacet->m_PVerts[0], *pFacet->m_PVerts[1]);
	double S2 = TriArea(p, *pFacet->m_PVerts[1], *pFacet->m_PVerts[2]);
	double S3 = TriArea(p, *pFacet->m_PVerts[2], *pFacet->m_PVerts[0]);
	double u = S1 / S; double v = S2 / S;  double w = S3 / S;
	VECTOR3D nor = vec[2] * u + vec[0] * v + vec[1] * w;
	nor.Normalize();
	return nor;
}
void  Ccalcubase::CalOnLineNorByFacNor(PHEDGE pHEdge, VECTOR3D& pNormal)
{
	double a, a1;
	PFACETTRI  pf = (PFACETTRI)pHEdge->pFacetAdj;
	PFACETTRI  pf1 = (PFACETTRI)pHEdge->pHEdgePair->pFacetAdj;
	double s = FacetTriArea(pf, a);
	double s1 = FacetTriArea(pf1, a1);
	double u = s / (s + s1);  double v = s1 / (s + s1);
	pNormal = (*pf->m_PFacetNorm) * u + (*pf1->m_PFacetNorm) * v;
	//pNormal = (*pf->m_PFacetNorm+ *pf1->m_PFacetNorm)/2.0;
	pNormal.Normalize();
}

bool Ccalcubase::IsHeight(POINT3D p0, POINT3D p1, POINT3D p2)
{
	VECTOR3D P0P2 = p2 - p0;
	VECTOR3D P1P0 = p0 - p1;
	VECTOR3D P1P2 = p2 - p1;
	double s = Ccalcubase::TriArea(p0, p1, p2);

	double dnext = P0P2.GetLength();
	double h = 2 * s / dnext;
	double angle = _AngleBetween3D(P1P0, P1P2);
	if ((h / dnext >= 0.3) || (angle < PI/2))
	{
		return true;
	}
	return false;
}