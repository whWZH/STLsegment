#include "StdAfx.h"
#include "SDF.h"
/////////////////////////////////////////////////////SDF聚类算法
double CSDF::SDF_calcu(PFACETTRI faceA,Octree* pOctree,vec_PFACETTRI& m_vecPFacetTri)
{  
	vector<VECTOR3D> m_vecVECTOR3D;
	vector<double> vec_SDF;
	POINT3D A=Ccalcubase::Center(faceA->pHEdge);
	VECTOR3D nA=(VECTOR3D)(*faceA->m_PFacetNorm);
	Ccalcubase::create_cone_SDF(m_vecVECTOR3D,A,nA,45,12);
	for (int ii=0;ii<m_vecVECTOR3D.size();ii++)
	{
		PFACETTRI getface=NULL;
		POINT3D cent;
		getface=pOctree->get_face_new(pOctree,m_vecVECTOR3D[ii],A,faceA,m_vecPFacetTri);
		if (getface!=NULL)
		{
			cent=Ccalcubase::Center(getface->pHEdge);
			double sdf=Ccalcubase::dis(A,cent);
			vec_SDF.push_back(sdf);
		}
	}
	Ccalcubase::SDF_delet(vec_SDF);
	double _SDF=0;
	for (int i=0;i<vec_SDF.size();i++)
	{
		_SDF=_SDF+vec_SDF[i];
	}
   return _SDF/vec_SDF.size();
}
void CSDF::SDF_kmeans(int k,vector<vector<PVERT>> &RBLOOP,vec_PFACETTRI& m_vecPFacetTri,vector<POINT3D>& vec_box,vector<vector<PFACETTRI>>& Rbox)
{
	Rbox.clear();
	map<double,PFACETTRI> map_SDF;
	vector<double> vec_SDF;
	vector<double> SDF_seg[10];
	vector<double> SDF_center;
	//set<double> set_sdf;
	set<double>::iterator SDF_it;
	///////////////8叉树
	Octree* pOctree;
	pOctree=NULL;
	//pOctree=new Octree;
	pOctree->creat_octree(pOctree,m_vecPFacetTri,vec_box);
	//////////////////////////////////计时
	time_t start,ends;
	start=clock();
	for (int i=0;i<m_vecPFacetTri.size();i++)
	{
		//double sdf=SDF_calcu(m_vecPFacetTri[i],pOctree,m_vecPFacetTri);
		double sdf=SDF_featurFACE(m_vecPFacetTri[i],pOctree,m_vecPFacetTri);
		sdf=sdf+(1-Ccalcubase::calcuPOINTV(m_vecPFacetTri[i]));
		map_SDF.insert(pair<double,PFACETTRI>(sdf,m_vecPFacetTri[i]));
		vec_SDF.push_back(sdf);
		//set_sdf.insert(sdf);
	}  
	ends=clock();
	double time1=difftime(ends,start);
	time1=0.001*time1;
	CString s;
	s.Format(_T("%f"),time1);
	AfxMessageBox(s);
	//////////////////////////////
	/*for (SDF_it=set_sdf.begin();SDF_it!=set_sdf.end();SDF_it++)
	{
		vec_SDF.push_back(*SDF_it);
	}*/
	sort(vec_SDF.begin(),vec_SDF.end());
	Ccalcubase::SDF_nomal(vec_SDF,map_SDF);
	/////////////////////////vec_SDF里面是排序之后的面表
	for (int i=1;i<=2*k;i++)
	{
		if ((i%2)!=0)
		{
			double num1=((double)i)/(2*k);
			double num2=num1*(vec_SDF.size()-1);
			SDF_center.push_back(vec_SDF[num2]);
		}
	}
	/////////////////////////k_means聚类
	while(1)
	{
		double num=0;
		//////////先分类
		for (int i=0;i<vec_SDF.size();i++)
		{
			int c=Ccalcubase::SDF_center(SDF_center,vec_SDF[i]);
			SDF_seg[c].push_back(vec_SDF[i]);
		}
		break;
		/////////每个类里找新聚点
		for (int i=0;i<k;i++)
		{
			double ki=0;
			ki=Ccalcubase::SDF_new_center(SDF_seg[i]);
			if (SDF_center[i]<=ki*1.1&&SDF_center[i]>=ki*0.9)
			{
				SDF_center[i]=ki;
				num++;
			}
			else
			{
				SDF_center[i]=ki;
			}
		}
		if (num==k)
		{
			break;
		}
	}
	for (int i=0;i<k;i++)
	{
		vector<PFACETTRI> vec_temp;
		for (int j=0;j<SDF_seg[i].size();j++)
		{
			vec_temp.push_back(map_SDF[SDF_seg[i][j]]);
		}
		Rbox.push_back(vec_temp);
	}
}
////////////////////////////////////////////////////////////////////////////SDF直方图
void CSDF::SDF_gram(int k,vector<vector<PVERT>> &RBLOOP,vector<double>& sdf_gram,vec_PFACETTRI& m_vecPFacetTri,vector<POINT3D>& vec_box,vec_PVERT& m_vecPVert)
{
	map<double,PFACETTRI> map_SDF;
	vec_PVERT sdf_point;
	vector<double> vec_SDF;
	//vector<int> sdf_gram;
	sdf_gram.resize(1000);
	double dis=0;
	/////////////////8叉树
	Octree* pOctree;
	pOctree=NULL;
	pOctree->creat_octree(pOctree,m_vecPFacetTri,vec_box);
	//////////////////////////////////计时
	time_t start,ends;
	start=clock();
	for (int i=0;i<m_vecPFacetTri.size();i++)
	{
		//double sdf=SDF_calcu(m_vecPFacetTri[i],pOctree,m_vecPFacetTri);
		double sdf=SDF_featurFACE(m_vecPFacetTri[i],pOctree,m_vecPFacetTri);
		//double sdf=(1-Ccalcubase::calcuPOINTV(m_vecPFacetTri[i]));
		//////////////////加上凹凸度
		sdf=sdf+(1-Ccalcubase::calcuPOINTV(m_vecPFacetTri[i]));
		//sdf=sdf*(1-Ccalcubase::calcuPOINTV(m_vecPFacetTri[i]));
		map_SDF.insert(pair<double,PFACETTRI>(sdf,m_vecPFacetTri[i]));
		vec_SDF.push_back(sdf);
	}  
	/*for (int i=0;i<m_vecPVert.size();i++)
	{
		double sdf=SDF_featurPOINT(m_vecPVert[i],pOctree,sdf_point,m_vecPFacetTri);
		sdf=sdf+(2-Ccalcubase::calcuNoV(m_vecPVert[i]));
		map_SDF.insert(pair<double,PVERT>(sdf,m_vecPVert[i]));
		vec_SDF.push_back(sdf);
	}*/
	sdf_point.clear();
	ends=clock();
	double time1=difftime(ends,start);
	time1=0.001*time1;
	CString s;
	s.Format(_T("%f"),time1);
	AfxMessageBox(s);
	////////////////////////////////
	sort(vec_SDF.begin(),vec_SDF.end());
    Ccalcubase::SDF_nomal(vec_SDF,map_SDF);
	//Ccalcubase::SDF_nomal_old(vec_SDF,map_SDF);
	///////////////////////////vec_SDF里面是排序之后的面表
	/////////////////////////////////////////////////SDF直方图
	dis=0.001;
	for (int i=0;i<vec_SDF.size();i++)
	{
		for (int j=1;j<=1000;j++)
		{
			if (vec_SDF[i]<=(j*dis)&&vec_SDF[i]>((j-1)*dis))
			{
				sdf_gram[j-1]=sdf_gram[j-1]+1;
			}
		}
	}
}
void CSDF::SDF_GMM(int k,vector<vector<PVERT>> &RBLOOP,vec_PFACETTRI& m_vecPFacetTri,vector<POINT3D>& vec_box,vector<vector<PFACETTRI>>& Rbox)
{
	map<double,PFACETTRI> map_SDF;
	vector<double> SDF_center;
	vector<double> vec_SDF;
	vector<double> SDF_seg[10];
	double dis=0;
	/////////////////8叉树
	Octree* pOctree;
	pOctree=NULL;
	pOctree->creat_octree(pOctree,m_vecPFacetTri,vec_box);
	//////////////////////////////////计时
	time_t start,ends;
	start=clock();
	for (int i=0;i<m_vecPFacetTri.size();i++)
	{
	    //double sdf=SDF_calcu(m_vecPFacetTri[i],pOctree,m_vecPFacetTri);
		///////加上凹凸度；
		double sdf=SDF_featurFACE(m_vecPFacetTri[i],pOctree,m_vecPFacetTri);
		//sdf=sdf+(1-Ccalcubase::calcuPOINTV(m_vecPFacetTri[i]));
		//sdf=sdf*(1-Ccalcubase::calcuPOINTV(m_vecPFacetTri[i]));
		map_SDF.insert(pair<double,PFACETTRI>(sdf,m_vecPFacetTri[i]));
		vec_SDF.push_back(sdf);
	}  
	ends=clock();
	double time1=difftime(ends,start);
	time1=0.001*time1;
	CString s;
	s.Format(_T("%f"),time1);
	AfxMessageBox(s);
	////////////////////////////////
	/////////////////////////////////////////析构八叉树
	///////////////////////////////////////////////////////////////////
	sort(vec_SDF.begin(),vec_SDF.end());
	Ccalcubase::SDF_nomal(vec_SDF,map_SDF);
	//Ccalcubase::SDF_nomal_old(vec_SDF,map_SDF);
	///////////////////////////vec_SDF里面是排序之后的面表
	for (int i=1;i<=2*k;i++)
	{
		if ((i%2)!=0)
		{
			double num1=((double)i)/(2*k);
			double num2=num1*(vec_SDF.size()-1);
			SDF_center.push_back(vec_SDF[num2]);
		}
	}
	/////////////////////////k_means聚类
	while(1)
	{
		double num=0;
		//////////先分类
		for (int i=0;i<vec_SDF.size();i++)
		{
			int c=Ccalcubase::SDF_center(SDF_center,vec_SDF[i]);
			SDF_seg[c].push_back(vec_SDF[i]);
		}
		/////////每个类里找新聚点
		for (int i=0;i<k;i++)
		{
			double ki=0;
			ki=Ccalcubase::SDF_new_center(SDF_seg[i]);
			if (SDF_center[i]<=ki*1.1&&SDF_center[i]>=ki*0.9)
			{
				SDF_center[i]=ki;
				num++;
			}
			else
			{
				SDF_center[i]=ki;
			}
		}
		if (num==k)
		{
			break;
		}
	}
	////////////////////////////////////////
	vector<double> vec_zmj;
	vector<double> vec_ct;
	for (int i=0;i<k;i++)
	{
		double temp=(double)(SDF_seg[i].size());
		double temp1=(double)(vec_SDF.size());
		double temp2=(double)temp/temp1;
		vec_zmj.push_back(temp2);
		vec_ct.push_back(1);
	}
	Ccalcubase::GMM_seg(vec_SDF,Rbox,map_SDF,k,SDF_center,vec_ct,vec_zmj);
}
/////////////////////////////////////////////////////////实验新算法
double CSDF::SDF_featurPOINT(PVERT& Apot,Octree* pOctree,vec_PVERT& sdf_point,vec_PFACETTRI& m_vecPFacetTri)
{
	double theSDF=0;
	POINT3D A=(POINT3D)(*Apot);
	VECTOR3D nA=Apot->Normal;
	nA.dx=-nA.dx;nA.dy=-nA.dy;nA.dz=-nA.dz;
	PFACETTRI getface=NULL;
	POINT3D cent;
	getface=pOctree->get_face_new1(pOctree,nA,A,(PFACETTRI)(Apot->pHEdgeOut->pFacetAdj),m_vecPFacetTri);
	if (getface!=NULL)
	{
		//PVERT pot=new VERT;
		VECTOR3D theDIS;
		cent=Ccalcubase::Center(getface->pHEdge);
		/*pot->x=(cent.x+A.x)*0.5;
		pot->y=(cent.y+A.y)*0.5;
		pot->z=(cent.z+A.z)*0.5;*/
		//sdf_point.push_back(pot);
		//pot->theEND=1;
		theDIS=cent-A;
		theSDF=theDIS.GetLength();
	}
	return theSDF;
}
void CSDF::SDF_getPOINT(vec_PFACETTRI& m_vecPFacetTri,vec_PVERT& m_vecPVert,vector<POINT3D>& vec_box)
{
	vec_PVERT sdf_point;
	/////////////////8叉树
	Octree* pOctree;
	pOctree=NULL;
	pOctree->creat_octree(pOctree,m_vecPFacetTri,vec_box);
	//////////////////////////////////计时
	time_t start,ends;
	start=clock();
	int nk=m_vecPVert.size();
	for (int i=0;i<nk;i++)
	{
		SDF_featurPOINT(m_vecPVert[i],pOctree,sdf_point,m_vecPFacetTri);
	}
	ends=clock();
	double time1=difftime(ends,start);
	time1=0.001*time1;
	CString s;
	s.Format(_T("%f"),time1);
	AfxMessageBox(s);
	delete pOctree;
	//pCalcubase->delet_octree(pOctree);
	//pCalcubase->segment_cloud(sdf_point,Rbox,m_vecPFacetTri);
}
double CSDF::SDF_featurFACE(PFACETTRI faceA,Octree* pOctree,vec_PFACETTRI& m_vecPFacetTri)
{
	double theSDF=0;
	POINT3D A=Ccalcubase::Center(faceA->pHEdge);
	VECTOR3D nA=(VECTOR3D)(*faceA->m_PFacetNorm);
	nA.dx=-nA.dx;nA.dy=-nA.dy;nA.dz=-nA.dz;
	PFACETTRI getface=NULL;
	POINT3D cent;
	vec_PFACETTRI vecface;
	Ccalcubase::FindPOneRFAC(faceA,vecface);
	double temSDF=0;
	for (int i=0;i<vecface.size();i++)
	{
			getface=pOctree->get_face_new1(pOctree,nA,A,faceA,m_vecPFacetTri);
	       if (getface!=NULL)
	    {
		   //PVERT pot=new VERT;
		  VECTOR3D theDIS;
		  cent=Ccalcubase::Center(getface->pHEdge);
		  /*pot->x=(cent.x+A.x)*0.5;
		  pot->y=(cent.y+A.y)*0.5;
		  pot->z=(cent.z+A.z)*0.5;*/
		  //sdf_point.push_back(pot);
		  //pot->theEND=1;
		  theDIS=cent-A;
		  theSDF=theDIS.GetLength();
	    }
		   temSDF=temSDF+theSDF;
	}
	return temSDF;
}