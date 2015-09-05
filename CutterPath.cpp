#include "stdafx.h"
#include "CutterPath.h"
CcutterPath::CcutterPath()
{

}

CcutterPath::~CcutterPath()
{

	int nNum, nCnt;
	nNum = m_vecPPolyPolygons.size();
	for (nCnt = 0; nCnt < nNum; nCnt++)
	{
		delete (m_vecPPolyPolygons[nCnt]);
	}
	m_vecPPolyPolygons.clear();
}
///////////////////////刀具轨迹
void CcutterPath::get_plane(vector<PVERT>& Bloop,double& A,double& B,double& C)
{
	MatrixXd P1(3,3),P2(3,3),P3(3,3);
	double x,y,z,xx,yy,zz,xy,yz,xz;
	x=y=z=xx=yy=zz=xy=yz=xz=0;
	for (int i=0;i<Bloop.size();i++)
	{
		x+=Bloop[i]->x;
		y+=Bloop[i]->y;
		z+=Bloop[i]->z;
		xx+=Bloop[i]->x*Bloop[i]->x;
		yy+=Bloop[i]->y*Bloop[i]->y;
		zz+=Bloop[i]->z*Bloop[i]->z;
		xy+=Bloop[i]->x*Bloop[i]->y;
		yz+=Bloop[i]->y*Bloop[i]->z;
		xz+=Bloop[i]->x*Bloop[i]->z;
	}
	P1(0,0)=xx;
	P1(0,1)=xy;
	P1(0,2)=xz;
	P1(1,0)=xy;
	P1(1,1)=yy;
	P1(1,2)=yz;
	P1(2,0)=xz;
	P1(2,1)=yz;
	P1(2,2)=zz;
	P2=P1.inverse();
	A=P2(0,0)*(-1)*(x)+P2(0,1)*(-1)*y+P2(0,2)*(-1)*z;
	B=P2(1,0)*(-1)*(x)+P2(1,1)*(-1)*y+P2(1,2)*(-1)*z;
	C=P2(2,0)*(-1)*(x)+P2(2,1)*(-1)*y+P2(2,2)*(-1)*z;
}
void CcutterPath::change(PFACETTRI& theFac)
{
	VECTOR3D V1, V2;
	V1 = *theFac->m_PVerts[0] - *theFac->m_PVerts[1];
	V2 = *theFac->m_PVerts[1] - *theFac->m_PVerts[2];
	VECTOR3D V3 = V1 * V2;
	V3.Normalize();
	double XZ = V3 | *theFac->m_PFacetNorm;
	if (XZ<0)
	{
		PVERT temPOT;
		temPOT = theFac->m_PVerts[0];
		theFac->m_PVerts[0] = theFac->m_PVerts[1];
		theFac->m_PVerts[1] = temPOT;
	}
}
void CcutterPath::cut_partilize(vector<vector<PFACETTRI>>& Rbox, vec_PFACETTRI& m_vecPFacetTri)
{
	for (int i = 0; i < m_vecPFacetTri.size(); i++)
	{
		m_vecPFacetTri[i]->bStatus = 0;
	}
	queue<PFACETTRI> Q;
	vec_PFACETTRI tempPath, vecFac;
	int numID = 0;
	while (1)
	{
		int num = 0;
		for (int i = 0; i < m_vecPFacetTri.size(); ++i)
		{
			if (m_vecPFacetTri[i]->bStatus == 0 && m_vecPFacetTri[i]->becut==0)
			{
				Q.push(m_vecPFacetTri[i]);
				num++;
				break;
			}
		}
		if (num == 0)
		{
			break;
		}
		while (!Q.empty())
		{
			Ccalcubase::FindPOneRFAC_NEW(Q.front(), vecFac);
			for (int i = 0; i < vecFac.size(); ++i)
			{
				if (vecFac[i]->bStatus == 0&&vecFac[i]->becut==0)
				{
					vecFac[i]->bStatus = 1;
					tempPath.push_back(vecFac[i]);
					Q.push(vecFac[i]);
				}
			}
			Q.front()->becut = 1;
			Q.pop();
		}
		for (int i = 0; i < tempPath.size(); i++)
		{
			tempPath[i]->ID = numID;
		}
		Rbox[numID]=tempPath;
		numID++;
		tempPath.clear();
	}
}
void CcutterPath::CUT_line(OBB& theBOX, PPOLYGON& firstPLN, vector<vec_VECTOR3D>& CCnomal, vec_PHEDGE &vec_PH, vector<vec_PFACETTRI>& RBOX)
{
	vec_VECTOR3D Cnomal;
	for (const PVERT &pot : *theBOX.pBloop)
	{
		Cnomal.push_back(pot->Normal);
		firstPLN->m_vecPnts.push_back(*pot);
		vec_PH.push_back(pot->pHEdgeOut);
	}
	CCnomal.push_back(Cnomal);
}
void CcutterPath::get_interPOT(PPOLYGON &firstPLN, PPOLYGON &nextPLN, vec_VECTOR3D &Cnomal, double &h, double &r, vec_PHEDGE &vec_PH, vec_PFACETTRI& facelist, VECTOR3D &Tdirection)
{
	//////VT为行距方向
	VECTOR3D V1,V2,Vij,Nij,VT;
	int numPLN=firstPLN->m_vecPnts.size();
	for (int i=0;i<numPLN;++i)
	{
		//////////////////构造偏置面
		vec_VECTOR3D interVECTOR;
		POINT3D PlanePOT;
		V1=firstPLN->m_vecPnts[(i+1)%numPLN]-firstPLN->m_vecPnts[i];
		if (i==0)
		{
			V2=firstPLN->m_vecPnts.back()-firstPLN->m_vecPnts[0];
		}
		else
		{
			V2=firstPLN->m_vecPnts[i-1]-firstPLN->m_vecPnts[i];
		}
		Vij=V1+V2;Vij.Normalize();Nij=Cnomal[i];
		VT=Nij*Vij;
		double bijiao = VT | Tdirection;
		if (bijiao>0)
		{
			VT.dx = -VT.dx; VT.dy = -VT.dy; VT.dz = -VT.dz;
		}
		PlanePOT = firstPLN->m_vecPnts[i] + Nij * 5;
	//	////////////////////计算行距
		double L=0;
		double p=0;
		POINT3D temp_POT=firstPLN->m_vecPnts[i];
		PFACETTRI TEMP_FAC=(PFACETTRI)(vec_PH[i]->pFacetAdj);
		if (firstPLN->m_bIsfirst==TRUE)
		{
		    p=Ccalcubase::CalPointOnLineCurv(VT,vec_PH[i],temp_POT);
		    p=1/p;
		}
		else
		{
			p=Ccalcubase::CalculatePoint_Curv(TEMP_FAC,temp_POT,VT);
			p=1/p;
		}
		if (p < 0)
		{
			L = sqrt((8 * r*h*p) / (p + r));
		}
		else
		{
			L = sqrt((8 * r*h*p) / (p - r));
		}
		/////////////////////////////////////移动平面上的点，并映射回曲面上
		PlanePOT=PlanePOT+VT*L*5;
		//nextPLN->m_vecPnts.push_back(PlanePOT);
		/*VECTOR3D invesNij;
		invesNij.dx = -Nij.dx; invesNij.dy = -Nij.dy; invesNij.dz = -Nij.dz;
		PlanePOT = PlanePOT + invesNij * 5;
		nextPLN->m_vecPnts.push_back(PlanePOT);*/
		///////////////////////////////////////////映射
		PFACETTRI jiaoFac = static_cast<PFACETTRI>(vec_PH[i]->pFacetAdj);
		vec_PFACETTRI tem_vecFac;
		//Ccalcubase::FindPOneRFAC(jiaoFac, tem_vecFac);
		Ccalcubase::FindPtwoRFAC(jiaoFac, tem_vecFac);
		for (int num = 0; num < tem_vecFac.size(); num++)
		{
			POINT3D JPOT;
			VECTOR3D invesNij;
			invesNij.dx = -Nij.dx; invesNij.dy = -Nij.dy; invesNij.dz = -Nij.dz;
			POINT3D facePOT[3];
			for (int k = 0; k < 3; k++)
			{
				facePOT[k] = *tem_vecFac[num]->m_PVerts[k];
			}
			if (Ccalcubase::IntersectTriangle(PlanePOT, invesNij, facePOT[0], facePOT[1], facePOT[2], JPOT)==TRUE)
			//if (Ccalcubase::IntersectLinePlaneYT(tem_vecFac[num], PlanePOT, invesNij, JPOT))
			{
				nextPLN->m_vecPnts.push_back(JPOT);
				facelist.push_back(tem_vecFac[num]);
				break;
			}
		}
		///////////////////////////////YT
		//POINT3D tempPOT;
		//OffsetPcc(facelist[i], Cnomal[i], firstPLN->m_vecPnts[i], tempPOT, h, r);
		//nextPLN->m_vecPnts.push_back(tempPOT);
	}	
}
void CcutterPath::get_CCpath(OBB& theBOX, vec_PPOLYPOLYGON&	m_vecPPolyPolygons, double h, double r, vector<vec_PFACETTRI>& RBOX)
{
	//////////////////////构造分割平面
	CPlane  cutPlane;
	vec_POINT3D theloop;
	for (const PVERT &pot : *theBOX.pBloop)
	{
		theloop.push_back(*pot);
	}
	double A, B, C;
	MatrixXd P1(3, 3), P2(3, 3), P3(3, 3);
	double x, y, z, xx, yy, zz, xy, yz, xz;
	x = y = z = xx = yy = zz = xy = yz = xz = 0;
	for (int i = 0; i < theloop.size(); i++)
	{
		x += theloop[i].x;
		y += theloop[i].y;
		z += theloop[i].z;
		xx += theloop[i].x*theloop[i].x;
		yy += theloop[i].y*theloop[i].y;
		zz += theloop[i].z*theloop[i].z;
		xy += theloop[i].x*theloop[i].y;
		yz += theloop[i].y*theloop[i].z;
		xz += theloop[i].x*theloop[i].z;
	}
	P1(0, 0) = xx;
	P1(0, 1) = xy;
	P1(0, 2) = xz;
	P1(1, 0) = xy;
	P1(1, 1) = yy;
	P1(1, 2) = yz;
	P1(2, 0) = xz;
	P1(2, 1) = yz;
	P1(2, 2) = zz;
	P2 = P1.inverse();
	A = P2(0, 0)*(-1)*(x)+P2(0, 1)*(-1)*y + P2(0, 2)*(-1)*z;
	B = P2(1, 0)*(-1)*(x)+P2(1, 1)*(-1)*y + P2(1, 2)*(-1)*z;
	C = P2(2, 0)*(-1)*(x)+P2(2, 1)*(-1)*y + P2(2, 2)*(-1)*z;
	VECTOR3D linor; linor.dx = A; linor.dy = B; linor.dz = C;
	cutPlane.a = A; cutPlane.b = B; cutPlane.c = C; cutPlane.d = 1;
	 cutPlane.m_Normal = linor;
	///////////////////////////////
	PPOLYPOLYGON  pPolyPolygon=new POLYPOLYGON;//整个轨迹
	vector<vec_VECTOR3D> CCnomal;//刀触点法矢
	PPOLYGON firstPLN=new POLYGON;//初始轨迹
	VECTOR3D Tdirection;//行距大方向
	update_direct(theBOX, Tdirection);

	vec_PHEDGE vec_PH;//轨迹所在半边
	CUT_line(theBOX,firstPLN,CCnomal,vec_PH,RBOX);//得到初始轨迹
	firstPLN->m_bIsfirst=1;
	pPolyPolygon->m_vecPPolygons.push_back(firstPLN);
	///////////////////获得下一条轨迹
	int num = 0;
	while (true)
	{
		PPOLYGON nextPLN = new POLYGON;
		vec_PFACETTRI nextFac;//轨迹所在面片
		get_interPOT(firstPLN, nextPLN, CCnomal.back(), h, r, vec_PH, nextFac, Tdirection);
		if (nextFac.size()==0)
		{
			break;
		}
		dele_interPOT(firstPLN, nextPLN, nextFac);
		vec_VECTOR3D cnomal;
		vec_PH.clear();
		ccpath_snooth(cutPlane, nextFac, nextPLN, cnomal,vec_PH);
		if (nextFac.size() == 0)
		{
			break;
		}
		CCnomal.push_back(cnomal);
		pPolyPolygon->m_vecPPolygons.push_back(nextPLN);
		firstPLN = nextPLN;
		++num;
	}
	m_vecPPolyPolygons.push_back(pPolyPolygon);
}
void CcutterPath::dele_interPOT(PPOLYGON firstPLN, PPOLYGON &nextPLN, vec_PFACETTRI& facelist)
{
	PPOLYGON tempPLN = new POLYGON;
	for (auto temp:firstPLN->m_vecPnts)
	{
		tempPLN->m_vecPnts.push_back(temp);
	}
	vector<POINT3D>::iterator Fit;
	vec_PFACETTRI::iterator FacIT;
	/////////////////去折线
	while (true)
	{
		int del = 0;
		for (int i = 1; i < nextPLN->m_vecPnts.size(); ++i)
		{
			VECTOR3D V1, V2;
			V1 = tempPLN->m_vecPnts[i] - tempPLN->m_vecPnts[i - 1];
			V2 = nextPLN->m_vecPnts[i] - nextPLN->m_vecPnts[i - 1];
			V1.Normalize(); V2.Normalize();
			double bijiao = (V1 | V2);
			if (bijiao <= 0)
			{
				Fit = tempPLN->m_vecPnts.begin() + i;
				tempPLN->m_vecPnts.erase(Fit);
				Fit = nextPLN->m_vecPnts.begin() + i;
				nextPLN->m_vecPnts.erase(Fit);
				FacIT = facelist.begin() + i;
				//facelist.erase(FacIT);
				--i;
				++del;
			}
		}
		if (del==0)
		{
			break;
		}
	}
	////////////////去波动
	while (true)
	{
		int del = 0;
		for (int i = 1; i < nextPLN->m_vecPnts.size(); ++i)  // 高度波动
		{
			if (Ccalcubase::IsHeight(nextPLN->m_vecPnts[i - 1], nextPLN->m_vecPnts[i], nextPLN->m_vecPnts[(i + 1) % nextPLN->m_vecPnts.size()]))
			{
				++del;
				Fit = nextPLN->m_vecPnts.begin() + i;
				nextPLN->m_vecPnts.erase(Fit);
				FacIT = facelist.begin() + i;
				//facelist.erase(FacIT);
				if (i < 2)
				{
					--i;
				}
				else
				{
					i = i - 2;
				}
			}
		}
		if (del==0)
			break;
	}
	
}
void CcutterPath::update_direct(OBB& theBOX, VECTOR3D &Tdirection)
{
	vector<PVERT> theLOOP = *theBOX.pBloop;
	POINT3D thePOT[3];
	for (int i = 0; i < theLOOP.size();++i)
	{
		thePOT[0].x = theLOOP[i]->x + thePOT[0].x;
		thePOT[0].y = theLOOP[i]->y + thePOT[0].y;
		thePOT[0].z = theLOOP[i]->z + thePOT[0].z;
	}
	for (int i = 0; i < 4;++i)
	{
		thePOT[1] = theBOX.theBox[i] + thePOT[1];
	}
	for (int i = 4; i < 8;++i)
	{
		thePOT[2] = theBOX.theBox[i] + thePOT[2];
	}
	thePOT[2] = thePOT[2] * 0.25;
	thePOT[1] = thePOT[1] * 0.25;
	thePOT[0] = thePOT[0] /theLOOP.size();
	VECTOR3D v1, v2;
	v1 = thePOT[0] - thePOT[1];
	v2 = thePOT[0] - thePOT[2];
	Tdirection = v1.GetLength()>v2.GetLength() ? v1 : v2;
}
void CcutterPath::ccpath_snooth(CPlane& cutPlane, vec_PFACETTRI& facelist, PPOLYGON &nextPLN, vec_VECTOR3D &Cnomal, vec_PHEDGE &vec_PH)
{
	Cnomal.clear();
	PPOLYGON tempPLN=new POLYGON;
	POINT3D center;
	for (auto &each :nextPLN->m_vecPnts )
	{
		tempPLN->m_vecPnts.push_back(each);
		center = center + each;
	}
	vec_PFACETTRI temFacelist;
	for (auto &each:facelist)
	{
		temFacelist.push_back(each);
	}
	/////////////////////////YT
	//for (int i = 0; i < tempPLN->m_vecPnts.size(); ++i)
	//{
	//	VECTOR3D V = tempPLN->m_vecPnts[(i + 1) % tempPLN->m_vecPnts.size()] - tempPLN->m_vecPnts[i];
	//	V.Normalize();
	//	VECTOR3D NV = V*(*temFacelist[i]->m_PFacetNorm);
	//	NV.Normalize();
	//	CPlane FP;////两轨迹点之间的平面0
	//	vec_POINT3D newPL;/////两轨迹点之间的点
	//	vec_PFACETTRI newFac;///////两轨迹点之间的面片
	//	//FP.m_Normal = NV; FP.m_PointIn = tempPLN->m_vecPnts[(i + 1) % tempPLN->m_vecPnts.size()]; FP.GetParameter();
	//	FP.BuildPlane(tempPLN->m_vecPnts[(i + 1) % tempPLN->m_vecPnts.size()], tempPLN->m_vecPnts[i], NV);
	//	FP.GetParameter();
	//	PHEDGE temPH[3]; POINT3D temPOT[3]; POINT3D INtemPOT;///////三个交点及其所在半边
	//	temPH[0] = temFacelist[i]->pHEdge; temPH[1] = temFacelist[i]->pHEdge->pHEdgeNext;
	//	temPH[2] = temFacelist[i]->pHEdge->pHEdgeNext->pHEdgeNext;
	//	PHEDGE firstPH ;
	//	for (int num = 0; num < 3;num++)
	//	{
	//		BOOL isin = FALSE;
	//		VECTOR3D V1;
	//		isin=Ccalcubase::IntersectLinePlane(temPH[num], FP, temPOT[num]);
	//		if (isin==TRUE)
	//		{
	//			V1 = temPOT[num] - tempPLN->m_vecPnts[i];
	//			V1.Normalize();
	//			double V1V = V1 | V;
	//			if (V1V>0)
	//			{
	//				newPL.push_back(temPOT[num]);
	//				firstPH = temPH[num];
	//				break;
	//			}
	//		}
	//	}
	//	////////////////继续循环，找到下面一系列交点
	//	if (!newPL.empty())
	//	{
	//		BOOL bInter;
	//		POINT3D InterP = newPL[0]; PHEDGE pHE = firstPH;
	//		PHEDGE nexPH = NULL; POINT3D PointPre = tempPLN->m_vecPnts[i];
	//		POINT3D InterP_Temp;
	//		do
	//		{
	//			bInter = FALSE;
	//			bInter = Ccalcubase::IntersectFacPlane(PointPre, pHE, InterP, FP, nexPH, InterP_Temp);
	//			if (bInter)
	//			{
	//				PointPre = InterP;
	//				pHE = nexPH;
	//				InterP = InterP_Temp;
	//				nexPH = NULL;
	//				VECTOR3D v3 = InterP_Temp - tempPLN->m_vecPnts[(i + 1) % tempPLN->m_vecPnts.size()];
	//				v3.Normalize();
	//				double v3v = V | v3;
	//				if (v3v >= 0)
	//				{
	//					break;
	//				}
	//				newPL.push_back(InterP_Temp);
	//				newFac.push_back(static_cast<PFACETTRI>(pHE->pFacetAdj));
	//			}
	//		} while (pHE->pFacetAdj != temFacelist[(i + 1) % temFacelist.size()]);
	//		nextPLN->m_vecPnts.insert(nextPLN->m_vecPnts.begin() + i, newPL.begin(), newPL.end());
	//		facelist.insert(facelist.begin() + i, newFac.begin(), newFac.end());
	//	}
	//}
	//////////////////////////////
	/////cut line
	/////////////////////////////////构造分割平面
	center.x = center.x / nextPLN->m_vecPnts.size(); center.y = center.y/ nextPLN->m_vecPnts.size();
	center.z = center.z/ nextPLN->m_vecPnts.size();
	cutPlane.m_PointIn = center;
	cutPlane.GetParameter();
	nextPLN->m_vecPnts.clear();
	facelist.clear();
	////////////////////////////找到首个交点
	PHEDGE firstPH;
	POINT3D PointPre;
	POINT3D InterP;
	PHEDGE pHE;
	for (int numi = 0; numi < temFacelist.size();numi++)
	{
		PHEDGE temPH[3]; POINT3D temPOT[3]; POINT3D INtemPOT;///////三个交点及其所在半边
		temPH[0] = temFacelist[numi]->pHEdge; temPH[1] = temFacelist[numi]->pHEdge->pHEdgeNext;
		temPH[2] = temFacelist[numi]->pHEdge->pHEdgeNext->pHEdgeNext;
		BOOL isbreak = FALSE;
		for (int num = 0; num < 3; num++)
		{
			BOOL isin = FALSE;
			VECTOR3D V1;
			isin = Ccalcubase::IntersectLinePlane(temPH[num], cutPlane, temPOT[num]);
			if (isin == TRUE)
			{
				pHE = temPH[num];
				InterP = temPOT[num];
				PointPre = tempPLN->m_vecPnts[0];
				nextPLN->m_vecPnts.push_back(InterP);
				facelist.push_back(static_cast<PFACETTRI>(pHE->pFacetAdj));
				vec_PH.push_back(pHE);
				VECTOR3D temv = Ccalcubase::CalcuPntInLineNormal(InterP, pHE);
				Cnomal.push_back(temv);
				isbreak = TRUE;
				break;
			}
		}
		if (isbreak==TRUE)
		{
			break;
		}
	}
	if (vec_PH.size()==0)
	{
		return ;
	}
	///////////////////////////
	BOOL bInter;
	PHEDGE nexPH; 
	POINT3D InterP_Temp;
	///////////////////////////////////////////////////////////////
	do
	{
		bInter = FALSE;
		bInter = Ccalcubase::IntersectFacPlane(PointPre, pHE, InterP, cutPlane, nexPH, InterP_Temp);
		if (bInter&&InterP_Temp != nextPLN->m_vecPnts[0])
		{
			PointPre = InterP;
			pHE = nexPH;
			InterP = InterP_Temp;
			nexPH = NULL;
			nextPLN->m_vecPnts.push_back(InterP_Temp);
			facelist.push_back(static_cast<PFACETTRI>(pHE->pFacetAdj));
			VECTOR3D temv=Ccalcubase::CalcuPntInLineNormal(InterP_Temp, pHE);
			Cnomal.push_back(temv);
			vec_PH.push_back(pHE);
		}
	} while (InterP_Temp != nextPLN->m_vecPnts[0]);
}

//PFACETTRI CcutterPath::OffsetPcc(PFACETTRI pFacet, VECTOR3D StepDir, POINT3D CurrentP, POINT3D& OffsetP, double ScallopHeight, double CutterR)
//{
//	CPlane plane;
//	VECTOR3D CurrentPnor;
//	VECTOR3D OffsetPnor;
//	//CalBase::CalculateNor(pFacet, CurrentP, CurrentPnor);
//	Ccalcubase::CalculatePointNor(pFacet, CurrentP, CurrentPnor);
//	double tempa = _AngleBetween3D(CurrentPnor, VECTOR3D(0, 1, 0));
//	double Interval = 0;
//	if (IS_ZERO(tempa - PI / 2.0))   //垂直方向
//	{
//		//double tempk = CalBase::CalculatePointCurv(pFacet, CurrentP, VECTOR3D(0, 1, 0));
//		double tempk = Ccalcubase::CalculatePoint_Curv(pFacet, CurrentP, VECTOR3D(0, 1, 0));
//		double kr = fabs(1 / tempk);
//		Interval = sqrt((8 * ScallopHeight * CutterR * kr) / (kr + CutterR));
//		OffsetP = CurrentP + VECTOR3D(0, 1, 0) * Interval;
//
//		// p所在面片
//		vec_PFACETTRI  NeighborFac;
//		//CalBase::FacetNeighborFacet(pstlmodel, pFacet, NeighborFac);
//		//CalBase::FacetNeighFacet(pFacet, NeighborFac);
//		Ccalcubase::Facet_nPNeighFacet(pFacet, NeighborFac);
//		for (unsigned int i = 0; i < NeighborFac.size(); i++)
//		{
//			if (Ccalcubase::PointInTri(NeighborFac.at(i), OffsetP))
//			{
//				PFACETTRI pF = NeighborFac.at(i);
//				Ccalcubase::CalculateNor(pF, OffsetP, OffsetPnor);
//				//vecNormal.push_back(OffsetPnor);
//				return pF;
//			}
//		}
//	}
//	else if (IS_ZERO(tempa))   // 水平方向
//	{
//		VECTOR3D vec0 = VECTOR3D(0, 1, 0) * StepDir;
//		//double tempk = CalBase::CalculatePointCurv(pFacet, CurrentP, vec0);
//		double tempk = Ccalcubase::CalculatePoint_Curv(pFacet, CurrentP, vec0);
//		double kr = fabs(1 / tempk);
//		Interval = sqrt((8 * ScallopHeight * CutterR * kr) / (kr + CutterR));
//		OffsetP = CurrentP + vec0 * Interval;
//		vec_PFACETTRI  NeighborFac;
//		//CalBase::FacetNeighborFacet(pstlmodel, pFacet, NeighborFac);
//		//CalBase::FacetNeighFacet(pFacet, NeighborFac);
//		Ccalcubase::Facet_nPNeighFacet(pFacet, NeighborFac);
//		for (unsigned int i = 0; i < NeighborFac.size(); i++)
//		{
//			if (Ccalcubase::PointInTri(NeighborFac.at(i), OffsetP))
//			{
//				PFACETTRI pF = NeighborFac.at(i);
//				Ccalcubase::CalculateNor(pF, OffsetP, OffsetPnor);
//				//vecNormal.push_back(OffsetPnor);
//				return pF;
//			}
//		}
//	}
//	//POINT3D ProjectP = ProjectionP(pFacet, StepDir, CurrentP, plane);   // 投影点
//	POINT3D ProjectP = ProjectionP_Sample(pFacet, StepDir, CurrentP, plane, ScallopHeight, CutterR, Interval);
//	// -----------------------------------
//	double testInterval = _DistOf(ProjectP, CurrentP);
//	VECTOR3D testv = ProjectP - CurrentP;
//	double testd = testv.GetLength();
//	// -----------------------------------
//	VECTOR3D LineDir = plane.m_Normal * (ProjectP - CurrentP);
//	LineDir.Normalize();
//
//	vec_PFACETTRI vecPFacet;
//	//CalBase::FacetNeighborFacet(pstlmodel, pFacet, vecPFacet);
//	//CalBase::FacetNeighFacet(pFacet, vecPFacet);
//	Ccalcubase::Facet_nPNeighFacet(pFacet, vecPFacet);
//	for (unsigned int i = 0; i < vecPFacet.size(); i++)
//	{
//		PFACETTRI pIntersectF = vecPFacet.at(i);
//		//bool Intersect = CalBase::LineTriIntersect(pIntersectF, ProjectP, LineDir, OffsetP);
//		//bool Intersect = CalBase::Line_TriIntersect(pIntersectF, ProjectP, LineDir, OffsetP);
//		bool Intersect = Ccalcubase::IntersectLinePlane(pIntersectF, ProjectP, LineDir, OffsetP);
//		if (Intersect)   // 直线与邻域面片相交
//		{
//			//OffsetPnor = CalBase::CalPointInTriNorByCore(pIntersectF, OffsetP);
//			Ccalcubase::CalculatePointNor(pIntersectF, OffsetP, OffsetPnor);
//			//vecNormal.push_back(OffsetPnor);
//			return pIntersectF;
//		}
//
//		//Intersect = CalBase::LineTriIntersect(pIntersectF, ProjectP, VECTOR3D(0, 0, 0)-LineDir, OffsetP);
//		//Intersect = CalBase::Line_TriIntersect(pIntersectF, ProjectP, VECTOR3D(0, 0, 0) - LineDir, OffsetP);
//		Intersect = Ccalcubase::IntersectLinePlane(pIntersectF, ProjectP, VECTOR3D(0, 0, 0) - LineDir, OffsetP);
//		if (Intersect)   // 直线与一阶邻域面片相交
//		{
//			//OffsetPnor = CalBase::CalPointInTriNorByCore(pIntersectF, OffsetP);
//			Ccalcubase::CalculatePointNor(pIntersectF, OffsetP, OffsetPnor);
//			//vecNormal.push_back(OffsetPnor);
//			return pIntersectF;
//		}
//	}
//
//	// 或逐点搜索
//	PHEDGE pHEdge;
//	PHEDGE pHE[3];
//	pHE[0] = pFacet->pHEdge;
//	pHE[1] = pHE[0]->pHEdgePair;
//	pHE[2] = pHE[1]->pHEdgeNext;
//	bool inter = false;
//	POINT3D p;
//	for (int k = 0; k < 3; k++)
//	{
//		inter = PlanewithEdge(pHE[k], plane, p);
//		if (inter)     // 平面与边相交
//		{
//			VECTOR3D v1 = p - CurrentP;    // CurrentP 当前刀触点   p 平面与边的交点
//			VECTOR3D v2 = ProjectP - CurrentP;    // ProjectP 投影点
//			if ((v1 | v2) > 0.0)      // 交点与投影点在同一侧
//			{
//				pHEdge = pHE[k];
//				break;
//			}
//		}
//	}
//
//	PFACETTRI pFac = (PFACETTRI)pHEdge->pHEdgePair->pFacetAdj;    // 错误，上一步，未求得半边
//	//inter = CalBase::LineTriIntersect(pFac, ProjectP, LineDir, OffsetP);
//	inter = Ccalcubase::Line_TriIntersect(pFac, ProjectP, LineDir, OffsetP);
//	//inter = CalBase::IntersectLinePlane(pFac, ProjectP, LineDir, OffsetP);
//	if (inter)   // 找到与直线相交的面片
//	{
//		//OffsetPnor = CalBase::CalPointInTriNorByCore(pFac, OffsetP);
//		Ccalcubase::CalculatePointNor(pFac, OffsetP, OffsetPnor);
//		//vecNormal.push_back(OffsetPnor);
//		return pFac;
//	}
//	else
//	{
//		do
//		{
//			// 找到该面片中与平面相交的边
//			pHEdge = pHEdge->pHEdgePair->pHEdgePre;
//			inter = PlanewithEdge(pHEdge, plane, p);
//			if (!inter)
//			{
//				pHEdge = pHEdge->pHEdgePair->pHEdgeNext;
//				inter = PlanewithEdge(pHEdge, plane, p);
//			}
//
//			pFac = (PFACETTRI)pHEdge->pHEdgePair->pFacetAdj;
//			//inter = CalBase::LineTriIntersect(pFac, ProjectP, LineDir, OffsetP);
//			//inter = CalBase::Line_TriIntersect(pFac, ProjectP, LineDir, OffsetP);
//			inter = Ccalcubase::IntersectLinePlane(pFac, ProjectP, LineDir, OffsetP);
//		} while (!inter);
//	}
//	//OffsetPnor = CalBase::CalPointInTriNorByCore(pFac, OffsetP);
//	Ccalcubase::CalculatePointNor(pFac, OffsetP, OffsetPnor);
//	//vecNormal.push_back(OffsetPnor);
//	return pFac;
//}
//
//POINT3D CcutterPath::ProjectionP_Sample(PFACETTRI pFacet, VECTOR3D StepDir, POINT3D CurrentP, CPlane& plane, double ScallopHeight, double CutterR, double Interval)
//{
//	//double CurrentPk;   // 当前点沿行距方向的曲率
//	VECTOR3D  CurrentPNor; // 
//	VECTOR3D IntervalDir;    // 行距方向
//	//PHEDGE pHEdge;	
//	double a = 0.0;   // 投影方向与行距方向的夹角
//	Ccalcubase::CalculatePointNor(pFacet, CurrentP, CurrentPNor);
//	double tempa = _AngleBetween3D(VECTOR3D(0, 1, 0), CurrentPNor);
//	if (tempa > PI / 2.0)
//	{
//		IntervalDir = VECTOR3D(0, 1, 0) * StepDir;
//		a = PI - tempa;
//	}
//	else   if (tempa < PI / 2.0)
//	{
//		IntervalDir = StepDir * VECTOR3D(0, 1, 0);
//		a = tempa;
//	}
//
//	// -----------------------------------------------------
//	// 求某一方向曲率时，由于InvervalDir在各顶点的切平面上的投影矢量不同，
//	// 故没办法进行简单的面积插值
//	VECTOR3D tp;
//	//CalBase::vecProject(CurrentPNor, IntervalDir, tp);
//	tp = IntervalDir - CurrentPNor * (IntervalDir | CurrentPNor);
//	double CurrentPk = Ccalcubase::CalculatePoint_Curv(pFacet, CurrentP, tp);   // 沿行距方向的曲率，近似
//
//	double kr = fabs(1.0 / CurrentPk);
//	//  此处为投影行距, 需判断曲面类型
//	if (CurrentPk < 0.0)  // 凸
//	{
//		Interval = cos(a) * sqrt((8 * ScallopHeight * CutterR * kr) / (kr + CutterR));
//	}
//	else
//	{
//		Interval = cos(a) * sqrt((8 * ScallopHeight * CutterR * kr) / (kr - CutterR));
//	}
//
//	POINT3D ProjectP = CurrentP + IntervalDir * Interval;     // 刀触点对应偏置点的对应投影点
//	VECTOR3D n = IntervalDir * VECTOR3D(0, 1, 0);    // 平面法矢：行距方向t与y的叉积
//	plane = plane.BuildPlane(n, CurrentP);
//	return ProjectP;
//}
//bool CcutterPath::PlanewithEdge(PHEDGE pHEdge, CPlane Plane, POINT3D& p)
//{
//	VECTOR3D  ab = *pHEdge->pVertEnd - *pHEdge->pHEdgePre->pVertEnd;
//	double t = (Plane.d - (Plane.m_Normal | (*pHEdge->pHEdgePre->pVertEnd - POINT3D(0, 0, 0)))) / (Plane.m_Normal| ab);
//	if (t >= 0.0 && t <= 1.0)
//	{
//		p = *pHEdge->pHEdgePre->pVertEnd + ab * t;
//		return true;
//	}
//	// t 趋于无穷，重合
//	return false;
//}