#include "stdafx.h"
#include "OBB.h"
using namespace Eigen;
using Eigen::MatrixXd;
void COBB::makeOBB(vec_PFACETTRI& seg_BOX,OBB &ob)
{
	MatrixXd allAT(3,3);
	double num=seg_BOX.size();
	VECTOR3D Uvector3d;
	Uvector3d.dx=0;Uvector3d.dy=0;Uvector3d.dz=0;

	//map<PVERT,VECTOR3D> nomaLIST;
	for (int i=0;i<seg_BOX.size();i++)
	{
		for (int j=0;j<3;j++)
		{
			ob.theCenter.x=ob.theCenter.x+seg_BOX[i]->m_PVerts[j]->x;
			ob.theCenter.y=ob.theCenter.y+seg_BOX[i]->m_PVerts[j]->y;
			ob.theCenter.z=ob.theCenter.z+seg_BOX[i]->m_PVerts[j]->z;

			/*VECTOR3D tempVECTOR;
			tempVECTOR=Ccalcubase::CalcuVerNormal(seg_BOX[i]->m_PVerts[j]);
			nomaLIST.insert(pair<PVERT,VECTOR3D>(seg_BOX[i]->m_PVerts[j],tempVECTOR));*/
		}
		Uvector3d.dx=Uvector3d.dx+seg_BOX[i]->m_PFacetNorm->dx;
		Uvector3d.dy=Uvector3d.dy+seg_BOX[i]->m_PFacetNorm->dy;
		Uvector3d.dz=Uvector3d.dz+seg_BOX[i]->m_PFacetNorm->dz;
	}
	
	ob.theCenter.x=ob.theCenter.x/(num*3);
	ob.theCenter.y=ob.theCenter.y/(num*3);
	ob.theCenter.z=ob.theCenter.z/(num*3);

	Uvector3d.dx=Uvector3d.dx/num;
	Uvector3d.dy=Uvector3d.dy/num;
	Uvector3d.dz=Uvector3d.dz/num;
	/////////////////////////////////求出协方差矩阵
	for (int i=0;i<seg_BOX.size();i++)
	{
		MatrixXd A(3,1),AT(1,3),ATA(3,3);
		A(0,0)=seg_BOX[i]->m_PFacetNorm->dx-Uvector3d.dx;
		A(1,0)=seg_BOX[i]->m_PFacetNorm->dy-Uvector3d.dy;
		A(2,0)=seg_BOX[i]->m_PFacetNorm->dz-Uvector3d.dz;
		AT=A.transpose();
		ATA=A*AT;
		allAT(0,0)=allAT(0,0)+ATA(0,0);allAT(0,1)=allAT(0,1)+ATA(0,1);
		allAT(0,2)=allAT(0,2)+ATA(0,2);
		allAT(1,0)=allAT(1,0)+ATA(1,0);allAT(1,1)=allAT(1,1)+ATA(1,1);
		allAT(1,2)=allAT(1,2)+ATA(1,2);
		allAT(2,0)=allAT(2,0)+ATA(2,0);allAT(2,1)=allAT(2,1)+ATA(2,1);
		allAT(2,2)=allAT(2,2)+ATA(2,2);
	}
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<3;j++)
		{
			allAT(i,j)=allAT(i,j)/num;
		}
	}
	////////求三个特征向量
	MatrixXcd A0,A1,A2;
	EigenSolver<MatrixXd> es(allAT);
	
	
	double valu[3],a(0),b(0),c(0);
	valu[0]=es.eigenvalues().col(0)(0).real();
	valu[1]=es.eigenvalues().col(0)(1).real();
	valu[2]=es.eigenvalues().col(0)(2).real();
	a=c=valu[0];
	for (int i=0;i<3;i++)
	{
		if (valu[i]>c)
		{
			c=valu[i];
		}
		if (valu[i]<a)
		{
			a=valu[i];
		}
	}
	for (int i=0;i<3;i++)
	{
		if (valu[i]==a)
		{
			A2=es.eigenvectors().col(i);
		}
		if (valu[i]!=a&&valu[i]!=c)
		{
			A1=es.eigenvectors().col(i);
		}
		if (valu[i]==c)
		{
			A0=es.eigenvectors().col(i);
		}
	}

	/*A0=es.eigenvectors().col(0);
	A1=es.eigenvectors().col(1);
	A2=es.eigenvectors().col(2);*/

	ob.theX.dx=A0(0).real();ob.theX.dy=A0(1).real();ob.theX.dz=A0(2).real();
	ob.theY.dx=A1(0).real();ob.theY.dy=A1(1).real();ob.theY.dz=A1(2).real();
	ob.theZ.dx=A2(0).real();ob.theZ.dy=A2(1).real();ob.theZ.dz=A2(2).real();
	//////////////////////////构建变化矩阵
	Ccalcubase::create_current_Frenet(ob.MA,ob.theZ,ob.theCenter);
	VECTOR3D initX;
	initX.dx=1;initX.dy=0;initX.dz=1;
	ob.theX.Normalize();ob.theY.Normalize();ob.theZ.Normalize();
	double ang=(initX|ob.theX);
	MATRIX3D turnX;
	Ccalcubase::create_current_matriX(turnX,ang);
	ob.MA=ob.MA*turnX;
	/////////////////////////////构建包围盒
	updateBOX(seg_BOX,ob);
	////////////////////////////标记分支
	set<int> IDsize;
	vec_PFACETTRI vecpFac;
	for (int i=0;i<seg_BOX.size();i++)
	{
		Ccalcubase::FindPOneRFAC_NEW(seg_BOX[i],vecpFac);
		for (int j=0;j<vecpFac.size();j++)
		{
			IDsize.insert(vecpFac[j]->ID);
		}
	}
	if (IDsize.size()==2)
	{
		ob.isBranch=TRUE;
	}
	else
	{
		ob.isBranch=FALSE;
	}
}
void COBB::updateBOX(vec_PFACETTRI  &seg_BOX,OBB &ob)
{
	vec_PFACETTRI new_BOX;
	for (int i=0;i<seg_BOX.size();i++)
	{
		PFACETTRI tempFAC=new FACETTRI;
		PVERT tempPOT[3];
		for (int j=0;j<3;j++)
		{
			tempPOT[j]=new VERT;
			double mm[4];
			mm[0]=seg_BOX[i]->m_PVerts[j]->x;mm[1]=seg_BOX[i]->m_PVerts[j]->y;mm[2]=seg_BOX[i]->m_PVerts[j]->z;
			mm[3]=1;
			tempPOT[j]->x=(ob.MA.A[0][0]*mm[0]+ob.MA.A[1][0]*mm[1]+ob.MA.A[2][0]*mm[2]+ob.MA.A[3][0]*mm[3]);
			tempPOT[j]->y=(ob.MA.A[0][1]*mm[0]+ob.MA.A[1][1]*mm[1]+ob.MA.A[2][1]*mm[2]+ob.MA.A[3][1]*mm[3]);
			tempPOT[j]->z=(ob.MA.A[0][2]*mm[0]+ob.MA.A[1][2]*mm[1]+ob.MA.A[2][2]*mm[2]+ob.MA.A[3][2]*mm[3]);
			tempFAC->m_PVerts[j]=tempPOT[j];
		}
		new_BOX.push_back(tempFAC);
	}
	int nNum = seg_BOX.size();
	if (nNum < 1)
	{
		return;
	}
	POINT3D POT[8];
	POT[0].x=new_BOX[1]->m_PVerts[0]->x;
	POT[0].y=new_BOX[1]->m_PVerts[0]->y;
	POT[0].z=new_BOX[1]->m_PVerts[0]->z;
	POT[7].x=new_BOX[1]->m_PVerts[0]->x;
	POT[7].y=new_BOX[1]->m_PVerts[0]->y;
	POT[7].z=new_BOX[1]->m_PVerts[0]->z;

	for (int n = 1; n< nNum; n++)
	{   
		for (int j = 1; j<3; j++){

			if ( new_BOX[n]->m_PVerts[j]->x> POT[7].x)
			{
				POT[7].x = new_BOX[n]->m_PVerts[j]->x;
			} 
			else
			{
				if (new_BOX[n]->m_PVerts[j]->x < POT[0].x)
				{
					POT[0].x= new_BOX[n]->m_PVerts[j]->x;
				}
			}

			if (new_BOX[n]->m_PVerts[j]->y >POT[7].y)
			{
				POT[7].y=new_BOX[n]->m_PVerts[j]->y ;
			} 
			else
			{
				if (new_BOX[n]->m_PVerts[j]->y < POT[0].y)
				{
					POT[0].y =new_BOX[n]->m_PVerts[j]->y ;
				}
			}

			if (new_BOX[n]->m_PVerts[j]->z >POT[7].z)
			{
				POT[7].z =new_BOX[n]->m_PVerts[j]->z;
			} 
			else
			{
				if (new_BOX[n]->m_PVerts[j]->z <POT[0].z)
				{
					POT[0].z =new_BOX[n]->m_PVerts[j]->z;
				}
			}
		}
	}
	////////////////////////////构造还原矩阵
	MatrixXd matriT(4,4),matriB(4,4);
	for (int i=0;i<4;i++)
	{
		for (int j=0;j<4;j++)
		{
			matriT(i,j)=ob.MA.A[i][j];
		}
	}
	matriB=matriT.inverse();
	////////////////////////////////////
	////父节点包围盒的八个顶点顺序
	//       4---7
	//      /|  /|
	//     5---6 |
	//     | 3-|-2
	//     |/  |/
	//     0---1
	//     */
	POT[1].x=POT[7].x;POT[1].y=POT[0].y;POT[1].z=POT[0].z;
	POT[2].x=POT[7].x;POT[2].y=POT[7].y;POT[2].z=POT[0].z;
	POT[3].x=POT[0].x;POT[3].y=POT[7].y;POT[3].z=POT[0].z;
	POT[4].x=POT[0].x;POT[4].y=POT[7].y;POT[4].z=POT[7].z;
	POT[5].x=POT[0].x;POT[5].y=POT[0].y;POT[5].z=POT[7].z;
	POT[6].x=POT[7].x;POT[6].y=POT[0].y;POT[6].z=POT[7].z;
	for (int i=0;i<8;i++)
	{
		double mm[4];
		mm[0]=POT[i].x;mm[1]=POT[i].y;mm[2]=POT[i].z;
		mm[3]=1;
		POT[i].x=(matriB(0,0)*mm[0]+matriB(1,0)*mm[1]+matriB(2,0)*mm[2]+matriB(3,0)*mm[3]);
		POT[i].y=(matriB(0,1)*mm[0]+matriB(1,1)*mm[1]+matriB(2,1)*mm[2]+matriB(3,1)*mm[3]);
		POT[i].z=(matriB(0,2)*mm[0]+matriB(1,2)*mm[1]+matriB(2,2)*mm[2]+matriB(3,2)*mm[3]);
	}

	for (int i=0;i<8;i++)
	{
		ob.theBox.push_back(POT[i]);
	}
}