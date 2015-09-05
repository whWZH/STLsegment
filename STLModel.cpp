#include "StdAfx.h"
#include "STLModel.h"

CSTLModel::CSTLModel(void)
{
	double m_avgGeod=0;
}


CSTLModel::~CSTLModel(void)
{
}
//struct CmpE   //map 的比较函数
//{             
//	bool operator()(const Edge& he1, const Edge& he2)
//	{if (he1.pVertStart< he2.pVertStart
//	|| (he1.pVertStart == he2.pVertStart && he1.pVertEnd < he2.pVertEnd))
//	return true;
//	else
//		return false;
//	}
//};
//map<VERT,PVERT>Vmap;
//pair<map<VERT,PVERT>::iterator,bool> Vmappoint;
//map<VERT,PVERT>::iterator iter; 
//map<EDGE,PEDGE,CmpE>mapedge;        
//pair<map<EDGE,PEDGE,CmpE>::iterator,bool> mapE;
//map<EDGE,PEDGE,CmpE>::iterator iterE; 

//BOOL CSTLModel::LoadSTLFile(LPCTSTR stlfile)
//{
//	FILE* file;
//	int type=0;
//	if((file=fopen(stlfile,"r"))==NULL)
//		return FALSE;
//
//	char str[80];
//	PHEDGE HEdg[3]={NULL,NULL,NULL};
//	PEDGE Edg[3]={NULL,NULL,NULL};
//	while (fscanf(file,"%s",str)==1)
//	{
//		if (strncmp(str,"normal",6)==0)
//		{
//			for (int i=0;i<3;i++)
//			{
//				HEdg[i]=new HEDGE;Edg[i]=new EDGE;
//			}
//			PFACETTRI ptri=new FACETTRI;
//			PVECTOR3D normal;
//			normal=new VECTOR3D;
//			fscanf(file,"%lf %lf %lf",&(normal->dx),&(normal->dy),&(normal->dz));
//			m_vecFacetNorm.push_back(*normal);
//			ptri->m_PFacetNorm=normal;
//
//			fscanf(file,"%*s %*s");
//			PVERT point=NULL;
//			for (int i=0; i<3; i++)
//			{
//				point=new VERT;
//				fscanf(file,"%*s %lf %lf %lf",&(point->x),&(point->y),&(point->z));
//				//bool used=false;
//				Vmappoint=Vmap.insert(pair<VERT,PVERT>(*point,point));
//				if (Vmappoint.second==0)
//				{ 
//					//used=true;
//					iter=Vmap.find(*point);
//					ptri->m_PVerts[i]=iter->second;
//					delete point;
//					point=NULL;
//				}
//				if (Vmappoint.second==1)
//				{
//					m_vecPVert.push_back(point);
//					ptri->m_PVerts[i]=m_vecPVert.back();
//					ptri->m_PVerts[i]->pHEdgeOut=HEdg[i];
//				}
//				HEdg[i]->pFacetAdj = ptri;
//				Edg[i]->pVertStart = ptri->m_PVerts[i];
//				Edg[i]->pHEdgeAdj = HEdg[i];
//
//			}
//			ptri->pHEdge=HEdg[0];
//			m_vecPFacetTri.push_back(ptri);	
//			HEdg[0]->pVertEnd = ptri->m_PVerts[1];
//			HEdg[1]->pVertEnd = ptri->m_PVerts[2];
//			HEdg[2]->pVertEnd = ptri->m_PVerts[0];
//
//			Edg[0]->pVertEnd = HEdg[0]->pVertEnd;
//			Edg[1]->pVertEnd = HEdg[1]->pVertEnd;
//			Edg[2]->pVertEnd = HEdg[2]->pVertEnd;
//
//			PEDGE temp=NULL;
//			for(int i=0;i<3;i++)
//			{   
//				temp=new EDGE;
//				temp->pVertStart=Edg[i]->pVertEnd;
//				temp->pVertEnd=Edg[i]->pVertStart;
//				temp->pHEdgeAdj=Edg[i]->pHEdgeAdj;
//				mapE=mapedge.insert(pair<EDGE,PEDGE>(*temp,temp));
//				if(mapE.second==0)
//				{   
//					iterE = mapedge.find(*temp); 
//					iterE->second->pHEdgeAdj ->pHEdgePair =Edg[i]->pHEdgeAdj;
//					Edg[i]->pHEdgeAdj->pHEdgePair=iterE->second->pHEdgeAdj;
//				}
//				else
//				{
//					mapE=mapedge.insert(pair<EDGE,PEDGE>(*Edg[i],Edg[i]));
//					m_vecPEdge.push_back(Edg[i]);
//				}
//			}
//
//			HEdg[0]->pHEdgeNext=HEdg[1];
//			HEdg[1]->pHEdgeNext=HEdg[2];
//			HEdg[2]->pHEdgeNext=HEdg[0];
//			HEdg[0]->pHEdgePre=HEdg[2];
//			HEdg[1]->pHEdgePre=HEdg[0];
//			HEdg[2]->pHEdgePre=HEdg[1];
//			m_vecPHEdge.push_back(HEdg[0]);
//			m_vecPHEdge.push_back(HEdg[1]);
//			m_vecPHEdge.push_back(HEdg[2]);
//			
//		}    
//	}
//	///////////////////////////////修补半边
//	for (int i=0;i<m_vecPHEdge.size();i++)
//	{
//		if (m_vecPHEdge[i]->pHEdgePair==NULL)
//		{
//			for (int j=0;j<m_vecPHEdge.size();j++)
//			{
//				if (m_vecPHEdge[j]->pVertEnd==m_vecPHEdge[i]->pHEdgePre->pVertEnd)
//				{
//					m_vecPHEdge[i]->pHEdgePair=m_vecPHEdge[j];
//					break;
//				}
//			}
//		}
//	}
//	///////////////////////////////
//	char title[80];
//	if (GetFileTitle(stlfile,title,80)==0)
//	{
//		SetName(title);
//	}
//	m_bModified = TRUE;
//	return TRUE;
//}
BOOL CSTLModel::LoadSTLFile(LPCTSTR stlfile)
{
	//CString strTime;         //测试程序运行时间
	//long t1 = GetTickCount();//程序运行前系统时间

	FILE* file = NULL;
	//int type=0;
	if((file = fopen(stlfile,"r")) == NULL)
		return FALSE;

	char str[80];
	                                 
	map<POINT3D,PVERT>mapPVert;                           //点映射，用于冗余点的去除
	pair<POINT3D,PVERT> pairP_PV;                         
	pair<map<POINT3D,PVERT>::iterator,bool> V_isRepeat;   //插入的返回值，用于判断点是否为重复点

	multimap<POINT3D,PHEDGE> mmapPHE;                //用于寻找伙伴半边
	pair<POINT3D,PHEDGE> pairP_PHE;
	multimap<POINT3D,PHEDGE>::iterator itVStartH,itVEndH;   //指向Key值为起点的pair,半边
	multimap<POINT3D,PHEDGE>::iterator itupperH;

	multimap<POINT3D,PEDGE> mmapPE;                //用于去除冗余边
	pair<POINT3D,PEDGE> pairP_PE;
	multimap<POINT3D,PEDGE>::iterator itVStart;   //指向Key值为起点的pair，边
	multimap<POINT3D,PEDGE>::iterator itupper;

	while(fscanf(file,"%s",str)==1)
	{
		if (strncmp(str,"normal",6)==0)
		{
			PFACETTRI pTri = NULL;
			pTri=new FACETTRI();
			PVECTOR3D pNor = NULL;
			pNor = new VECTOR3D();

			fscanf(file,"%lf %lf %lf",&(pNor->dx),&(pNor->dy),&(pNor->dz));
            pTri->m_PFacetNorm = pNor;
			//m_vecFacetNorm.push_back(pNor);

			fscanf(file,"%*s %*s");

			for (int i=0;i<3;i++)
			{
				PVERT  pVer = NULL;
				pVer=new VERT();
				fscanf(file,"%*s %lf %lf %lf",&(pVer->x),&(pVer->y),&(pVer->z));
                
				POINT3D Pt = *pVer;
				pairP_PV.first = Pt;
				pairP_PV.second = pVer;
				V_isRepeat = mapPVert.insert(pairP_PV);
				if (V_isRepeat.second)                           //判断插入是否成功：1为插入成功，不是重复点；
				{                                                                //：0为插入不成功，是重复点
                    pTri->m_PVerts[i] = pVer;
					m_vecPVert.push_back(pVer);  //*****点表
				} 
				else
				{
                   pTri->m_PVerts[i] = V_isRepeat.first->second;                //将当前迭代器所指点的值给Vec
					delete pVer;                                 //释放动态空间，防止内存泄露 
					pVer = 0;
				}
			}
       //**************************建立拓扑关系****************************//
       //*****建立半边表
         PHEDGE pHE_AdjF[3];                                     //定义一个存放每个三角面片三条半边的数组指针
		 pHE_AdjF[0] = NULL;
		 pHE_AdjF[1] = NULL;
		 pHE_AdjF[2] = NULL;
	     pHE_AdjF[0] = new HEDGE();
		 pHE_AdjF[1] = new HEDGE();
		 pHE_AdjF[2] = new HEDGE();
		 PHEDGE pHE = NULL;

		 for (int j=0;j<3;j++)           
		 {
           pHE = pHE_AdjF[j];
           pHE->pVertEnd = pTri->m_PVerts[(j+1)%3];                //*半边的终点,0号半边的终点是1号
		   pHE->pHEdgeNext = pHE_AdjF[(j+1)%3];                    //*半边沿面片的下个半边
		   pHE->pHEdgePre = pHE_AdjF[(j+2)%3];                     //*半边沿面片的上一个半边
		   pHE->pFacetAdj = pTri;                                  //*半边相邻的面片

           pHE->pVertEnd->pHEdgeOut = pHE->pHEdgeNext;    //*半边j的终点向外发散的一个半边为半边j的下一条半边          
		   //*寻找伙伴半边 
		   //查找方法：找到以原半边的终点为起点的所有半边，然后在这些半边中查找半边终点等于原半边起点的半边
		   PHEDGE pHEdge = NULL;
		   PVERT pStart = NULL,pEnd = NULL;
		   pHEdge = pHE;
		   pStart = pTri->m_PVerts[j];
		   pEnd   = pHEdge->pVertEnd;

		   pairP_PHE.first = *pStart;                      //以半边的起点作为Key值
		   pairP_PHE.second = pHEdge;

		   if (mmapPHE.empty() != true)                    
		   {
			   itVStartH = mmapPHE.lower_bound(*pEnd);       //指向Key值以pHEdge终点作为起点的第一个pair
			   if (itVStartH != mmapPHE.end())
			   {
			   	   itupperH = mmapPHE.upper_bound(*pEnd);     //指向Key值以pHEdge终点作为起点的最后一个pair的后一个
				   while(itVStartH != itupperH)
				   {
					   if(itVStartH->second->pHEdgePair ==NULL)//如果该半边还没找到伙伴半边
					   {
						   if(itVStartH->second->pVertEnd == pStart)  //找到终点等于pHEdge起点的半边为其伙伴半边
						   {
							   pHEdge->pHEdgePair = itVStartH->second;
							   itVStartH->second->pHEdgePair = pHEdge;  //两半边互为伙伴半边
							   mmapPHE.erase(itVStartH);                  //删除已找到伙伴半边的map,2013.6.25添
							   break;
						   }
					   }
					   itVStartH++;
				   }
			   }
		   }
		   m_vecPHEdge.push_back(pHE);             //*****半边表  
		   if (pHEdge->pHEdgePair == NULL)   //假如没找到伙伴半边，就把半边添加到map以等待找伙伴半边，2013.6.25改            
		   {
			   mmapPHE.insert(pairP_PHE);
		   //*****建立边表
			   PEDGE pEdge = NULL;
			   pEdge = new EDGE();

			   PVERT pPreHEndV,pHEndV;
			   pPreHEndV = pTri->m_PVerts[j];         //半边i 的前条半边的终点是 半边i 的起点
			   pHEndV    = pHE->pVertEnd;
			   if (*pPreHEndV < *pHEndV)                                 //比较半边的起点和终点，
			   {   pEdge->pVertStart = pPreHEndV;                //*将较小的点作为边的起点，
				   pEdge->pVertEnd   = pHEndV;                   //*较大的作为边的终点
				} 
				else
				{	pEdge->pVertStart = pHEndV;
				    pEdge->pVertEnd   = pPreHEndV;	
				 }
				 pEdge->pHEdgeAdj  = pHE;            //*将当前的半边作为边相邻的一个半边

				   //pairP_PE.first = *pEdge->pVertStart; //以边的起点作为Key值
				   //pairP_PE.second = pEdge;

				   //   if (mmapPE.empty() != true)
				   //   {
				   //   itVStart = mmapPE.lower_bound(pairP_PE.first);     //指向Key值为起点的第一个pair
				   //   if (itVStart != mmapPE.end())
				   //   {
				   //      itupper = mmapPE.upper_bound(pairP_PE.first);   //指向Key值为起点的最后一个pair的后一个
				   //         while(itVStart != itupper)
				   //         {
				   //          if(((*itVStart).second->pVertEnd) == (pEdge->pVertEnd))
				   //          {
				   //	             delete pEdge;
				   //           pEdge=0;
				   //           break;
				   //          }
				   //          itVStart++;
				   //         }
				   //           }
				   //       }
				   //mmapPE.insert(pairP_PE);
				   /*	        if (pEdge!=0)*/
				   m_vecPEdge.push_back(pEdge);       //*****边表	
		   }
		 }
		 pTri->pHEdge = pHE;                                            //*将三角面片的最后一个半边给面片	
		 m_vecPFacetTri.push_back(pTri);   //*****面表			
		} 
	}

	fclose(file);    //2013.8.15,需要将文件关闭，一定要注意，否则会造成文件所占用内存泄露和在下次访问文件时出现问题
	file = NULL;     //2013.8.15,需要将文件指针指向空，这样做会防止出现游离指针，而对整个工程造成不必要的麻烦
	for (int i=0;i<m_vecPVert.size();i++)
	{
		m_vecPVert[i]->Normal=Ccalcubase::CalcuVerNormal(m_vecPVert[i]);
		m_vecPVert[i]->ID=i;
	}
	//itVStartH = mmapPHE.begin(); //2013.7.11	
	//while(itVStartH != mmapPHE.end())
	//{		
	//	POINT3D  PointEnd;
	//	BOOL isFind = FALSE;
	//	PointEnd = (POINT3D)*(itVStartH->second->pVertEnd);
	//	itVEndH = mmapPHE.find(PointEnd);

	//	if (itVEndH != mmapPHE.end())
	//	{
	//		while(itVEndH != mmapPHE.upper_bound(PointEnd))
	//		{
	//			if ((POINT3D)*(itVEndH->second->pVertEnd) == (POINT3D)(itVStartH->first))
	//			{
	//				itVStartH->second->pHEdgePair = itVEndH->second;
	//				itVEndH->second->pHEdgePair = itVStartH->second;
	//				mmapPHE.erase(itVStartH);
	//				mmapPHE.erase(itVEndH);
	//				itVStartH = mmapPHE.begin();
	//				isFind = TRUE;
	//				break;
	//			}
	//			itVEndH++;
	//		}
	//	}
	//	if (isFind == FALSE)	itVStartH++;
	//}

	//long t2 = GetTickCount();         //程序运行后系统时间
	//strTime.Format("time:%dms",t2-t1);//前后时间差为程序运行时间
	//AfxMessageBox(strTime);           //用消息框显示
	m_bModified = TRUE;
	return TRUE;
}
void CSTLModel::UpdateBox()
{
	int nNum = m_vecPFacetTri.size();
	if (nNum < 1)
	{
		return;
	}
	BOX3D one;
	one.x0=m_vecPFacetTri[1]->m_PVerts[0]->x;
	one.y0=m_vecPFacetTri[1]->m_PVerts[0]->y;
	one.z0=m_vecPFacetTri[1]->m_PVerts[0]->z;
	one.x1=m_vecPFacetTri[1]->m_PVerts[0]->x;
	one.y1=m_vecPFacetTri[1]->m_PVerts[0]->y;
	one.z1=m_vecPFacetTri[1]->m_PVerts[0]->z;
	
	for (int n = 1; n< nNum; n++)
	{   
		for (int j = 1; j<3; j++){

			if ( m_vecPFacetTri[n]->m_PVerts[j]->x> one.x1)
			{
				one.x1 = m_vecPFacetTri[n]->m_PVerts[j]->x;
			} 
			else
			{
				if (m_vecPFacetTri[n]->m_PVerts[j]->x < one.x0)
				{
					one.x0= m_vecPFacetTri[n]->m_PVerts[j]->x;
				}
			}

			if (m_vecPFacetTri[n]->m_PVerts[j]->y >one.y1)
			{
				one.y1=m_vecPFacetTri[n]->m_PVerts[j]->y ;
			} 
			else
			{
				if (m_vecPFacetTri[n]->m_PVerts[j]->y < one.y0)
				{
					one.y0 =m_vecPFacetTri[n]->m_PVerts[j]->y ;
				}
			}

			if (m_vecPFacetTri[n]->m_PVerts[j]->z >one.z1)
			{
				one.z1 =m_vecPFacetTri[n]->m_PVerts[j]->z;
			} 
			else
			{
				if (m_vecPFacetTri[n]->m_PVerts[j]->z <one.z0)
				{
					one.z0 =m_vecPFacetTri[n]->m_PVerts[j]->z;
				}
			}
		}
	}
	m_Box=one;
	POINT3D theMAX,theMIN;
	theMAX.x=one.x1;theMAX.y=one.y1;theMAX.z=one.z1;
	theMIN.x=one.x0;theMIN.y=one.y0;theMIN.z=one.z0;
	vec_box.clear();
	vec_box.push_back(theMIN);
	vec_box.push_back(theMAX);
}

void CSTLModel::Draw(COpenGLDC* pDC)
{
	pDC->SetMaterialColor(RGB(225,175,22));
	pDC->SetColor(RGB(255,255,0));
	int nNumF = m_vecPFacetTri.size();
	for(int i=0;i<nNumF;i++)
	{
		pDC->SetMaterialColor(RGB(225,175,22));
		pDC->DrawTriChip(m_vecPFacetTri[i]->m_PFacetNorm->dx,m_vecPFacetTri[i]->m_PFacetNorm->dy,m_vecPFacetTri[i]->m_PFacetNorm->dz,
			m_vecPFacetTri[i]->m_PVerts[0]->x,m_vecPFacetTri[i]->m_PVerts[0]->y,m_vecPFacetTri[i]->m_PVerts[0]->z,
			m_vecPFacetTri[i]->m_PVerts[1]->x,m_vecPFacetTri[i]->m_PVerts[1]->y,m_vecPFacetTri[i]->m_PVerts[1]->z,
			m_vecPFacetTri[i]->m_PVerts[2]->x,m_vecPFacetTri[i]->m_PVerts[2]->y,m_vecPFacetTri[i]->m_PVerts[2]->z);
	}
}

void CSTLModel::Draw2(COpenGLDC* pDC)
{
	COLORREF clr,clrold;
	clr = RGB(0,0,0);
	pDC->GetMaterialColor(clrold);
	pDC->SetMaterialColor(clr);
	
	for (int i=0;i<m_vecPVert.size();i++)
	{
		if (m_vecPVert[i]->bused==1)
		{
			pDC->SetColor(RGB(0,0,0));
			pDC->DrawPoint(*m_vecPVert[i]);
		}
		/*if (m_vecPVert[i]->theEND==1)
		{
			pDC->SetColor(RGB(255,255,255));
			pDC->DrawPoint(*m_vecPVert[i]);
		}*/
	}
	for (int i=0;i<m_vecPHEdge.size();i++)
	{
		if (m_vecPHEdge[i]->bused==1)
		{
			pDC->SetColor(RGB(0,0,0));
			pDC->DrawLine(*(m_vecPHEdge[i]->pHEdgePair->pVertEnd),*(m_vecPHEdge[i]->pVertEnd));
		}
	}
	pDC->SetMaterialColor(clrold);
}

bool CSTLModel::IsEmpty()
{
	return m_vecPFacetTri.empty();
}
void CSTLModel::GetBox(double x0,double y0,double z0,double x1,double y1,double z1)
{
	
	UpdateBox();
	m_one.x0=m_Box.x0;
	m_one.y0=m_Box.y0;
	m_one.z0=m_Box.z0;
	m_one.x1=m_Box.x1;
	m_one.y1=m_Box.y1;
	m_one.z1=m_Box.z1;
}

void CSTLModel::Draw1(COpenGLDC* pDC,vector<vector<PFACETTRI>>& Rbox)
{
	COLORREF old_clr; 
	for (int i=0;i<Rbox.size();i++)
	{
		pDC->GetMaterialColor(old_clr);
		pDC->SetMaterialColor(RGB((47*(i+1))%255,(183*i)%255,(197*i)%255));
		pDC->SetColor(RGB((47*(i+1))%255,(183*i)%255,(197*i)%255));
		for (int j=0;j<Rbox[i].size();j++)
		{
			pDC->DrawTriChip(Rbox[i][j]->m_PFacetNorm->dx,Rbox[i][j]->m_PFacetNorm->dy,Rbox[i][j]->m_PFacetNorm->dz,
				Rbox[i][j]->m_PVerts[0]->x,Rbox[i][j]->m_PVerts[0]->y,Rbox[i][j]->m_PVerts[0]->z,
				Rbox[i][j]->m_PVerts[1]->x,Rbox[i][j]->m_PVerts[1]->y,Rbox[i][j]->m_PVerts[1]->z,
				Rbox[i][j]->m_PVerts[2]->x,Rbox[i][j]->m_PVerts[2]->y,Rbox[i][j]->m_PVerts[2]->z);
		}
	}
}
void CSTLModel::Draw3(COpenGLDC* pDC,vector<vector<PFACETTRI>>& Rbox)
{
	COLORREF old_clr; 
	for (int i=0;i<Rbox.size();i++)
	{
		pDC->GetMaterialColor(old_clr);
		pDC->SetMaterialColor(RGB(200-20*i,255-15*i,255-22*i));
		pDC->SetColor(RGB(200-20*i,255-15*i,255-22*i));
		for (int j=0;j<Rbox[i].size();j++)
		{
			pDC->DrawTriChip(Rbox[i][j]->m_PFacetNorm->dx,Rbox[i][j]->m_PFacetNorm->dy,Rbox[i][j]->m_PFacetNorm->dz,
				Rbox[i][j]->m_PVerts[0]->x,Rbox[i][j]->m_PVerts[0]->y,Rbox[i][j]->m_PVerts[0]->z,
				Rbox[i][j]->m_PVerts[1]->x,Rbox[i][j]->m_PVerts[1]->y,Rbox[i][j]->m_PVerts[1]->z,
				Rbox[i][j]->m_PVerts[2]->x,Rbox[i][j]->m_PVerts[2]->y,Rbox[i][j]->m_PVerts[2]->z);
		}
	}
}
void CSTLModel::Draw_obb(COpenGLDC* pDC,vector<OBB>& vec_obb)
{
	for (int i=0;i<vec_obb.size();i++)
	{
		POINT3D xPt,yPt,zPt;
		xPt=vec_obb[i].theCenter+vec_obb[i].theX;
		yPt=vec_obb[i].theCenter+vec_obb[i].theY;
		zPt=vec_obb[i].theCenter+vec_obb[i].theZ;
		COLORREF old_clr;
		pDC->GetColor(old_clr);

		/*pDC->SetColor(RGB(255,0,0));
		pDC->DrawLine(vec_obb[i].theCenter,xPt);

		pDC->SetColor(RGB(0,255,0));
		pDC->DrawLine(vec_obb[i].theCenter,yPt);*/

		/*pDC->SetColor(RGB(0,0,255));
		pDC->DrawLine(vec_obb[i].theCenter,zPt);*/

		pDC->SetColor(RGB(0,0,0));
		pDC->DrawLine(vec_obb[i].theBox[0],vec_obb[i].theBox[1]);
		pDC->DrawLine(vec_obb[i].theBox[1],vec_obb[i].theBox[2]);
		pDC->DrawLine(vec_obb[i].theBox[2],vec_obb[i].theBox[3]);
		pDC->DrawLine(vec_obb[i].theBox[3],vec_obb[i].theBox[0]);
		pDC->DrawLine(vec_obb[i].theBox[5],vec_obb[i].theBox[6]);
		pDC->DrawLine(vec_obb[i].theBox[6],vec_obb[i].theBox[7]);
		pDC->DrawLine(vec_obb[i].theBox[7],vec_obb[i].theBox[4]);
		pDC->DrawLine(vec_obb[i].theBox[4],vec_obb[i].theBox[5]);
		pDC->DrawLine(vec_obb[i].theBox[0],vec_obb[i].theBox[5]);
		pDC->DrawLine(vec_obb[i].theBox[1],vec_obb[i].theBox[6]);
		pDC->DrawLine(vec_obb[i].theBox[2],vec_obb[i].theBox[7]);
		pDC->DrawLine(vec_obb[i].theBox[3],vec_obb[i].theBox[4]);

	    pDC->SetColor(old_clr); 
	}
	
}
void CSTLModel::Draw_ccPath(COpenGLDC* pDC,vec_PPOLYPOLYGON&	m_vecPPolyPolygons)
{
	POINT3D StarP,EndP;
	COLORREF clr,clrold;
	clr = RGB(0,0,0);
	pDC->GetMaterialColor(clrold);
	pDC->SetMaterialColor(clr);
	PPOLYGON pPoly;
	for (int i=0;i<m_vecPPolyPolygons.size();i++)
	{
		for (int j=0;j<m_vecPPolyPolygons[i]->m_vecPPolygons.size();j++)
		{
			pPoly = m_vecPPolyPolygons[i]->m_vecPPolygons[j];
			for (int k=0;k<pPoly->m_vecPnts.size()-1;k++)
			{
				StarP = pPoly->m_vecPnts[k];
				EndP  = pPoly->m_vecPnts[k+1];
				pDC->DrawLine(StarP,EndP);
			}
			StarP = pPoly->m_vecPnts[pPoly->m_vecPnts.size()-1];
			EndP  = pPoly->m_vecPnts[0];
			pDC->DrawLine(StarP,EndP);
		}
		/////////////test
		/*for (int j = 0; j < m_vecPPolyPolygons[i]->m_vecPPolygons[0]->m_vecPnts.size(); ++j)
		{
			StarP = m_vecPPolyPolygons[i]->m_vecPPolygons[0]->m_vecPnts[j];
			EndP = m_vecPPolyPolygons[i]->m_vecPPolygons[1]->m_vecPnts[j];
			pDC->DrawLine(StarP, EndP);
		}*/
	}
}