#include "StdAfx.h"
#include "STLModel.h"

CSTLModel::CSTLModel(void)
{
	double m_avgGeod=0;
}


CSTLModel::~CSTLModel(void)
{
}
//struct CmpE   //map �ıȽϺ���
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
//	///////////////////////////////�޲����
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
	//CString strTime;         //���Գ�������ʱ��
	//long t1 = GetTickCount();//��������ǰϵͳʱ��

	FILE* file = NULL;
	//int type=0;
	if((file = fopen(stlfile,"r")) == NULL)
		return FALSE;

	char str[80];
	                                 
	map<POINT3D,PVERT>mapPVert;                           //��ӳ�䣬����������ȥ��
	pair<POINT3D,PVERT> pairP_PV;                         
	pair<map<POINT3D,PVERT>::iterator,bool> V_isRepeat;   //����ķ���ֵ�������жϵ��Ƿ�Ϊ�ظ���

	multimap<POINT3D,PHEDGE> mmapPHE;                //����Ѱ�һ����
	pair<POINT3D,PHEDGE> pairP_PHE;
	multimap<POINT3D,PHEDGE>::iterator itVStartH,itVEndH;   //ָ��KeyֵΪ����pair,���
	multimap<POINT3D,PHEDGE>::iterator itupperH;

	multimap<POINT3D,PEDGE> mmapPE;                //����ȥ�������
	pair<POINT3D,PEDGE> pairP_PE;
	multimap<POINT3D,PEDGE>::iterator itVStart;   //ָ��KeyֵΪ����pair����
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
				if (V_isRepeat.second)                           //�жϲ����Ƿ�ɹ���1Ϊ����ɹ��������ظ��㣻
				{                                                                //��0Ϊ���벻�ɹ������ظ���
                    pTri->m_PVerts[i] = pVer;
					m_vecPVert.push_back(pVer);  //*****���
				} 
				else
				{
                   pTri->m_PVerts[i] = V_isRepeat.first->second;                //����ǰ��������ָ���ֵ��Vec
					delete pVer;                                 //�ͷŶ�̬�ռ䣬��ֹ�ڴ�й¶ 
					pVer = 0;
				}
			}
       //**************************�������˹�ϵ****************************//
       //*****������߱�
         PHEDGE pHE_AdjF[3];                                     //����һ�����ÿ��������Ƭ������ߵ�����ָ��
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
           pHE->pVertEnd = pTri->m_PVerts[(j+1)%3];                //*��ߵ��յ�,0�Ű�ߵ��յ���1��
		   pHE->pHEdgeNext = pHE_AdjF[(j+1)%3];                    //*�������Ƭ���¸����
		   pHE->pHEdgePre = pHE_AdjF[(j+2)%3];                     //*�������Ƭ����һ�����
		   pHE->pFacetAdj = pTri;                                  //*������ڵ���Ƭ

           pHE->pVertEnd->pHEdgeOut = pHE->pHEdgeNext;    //*���j���յ����ⷢɢ��һ�����Ϊ���j����һ�����          
		   //*Ѱ�һ���� 
		   //���ҷ������ҵ���ԭ��ߵ��յ�Ϊ�������а�ߣ�Ȼ������Щ����в��Ұ���յ����ԭ������İ��
		   PHEDGE pHEdge = NULL;
		   PVERT pStart = NULL,pEnd = NULL;
		   pHEdge = pHE;
		   pStart = pTri->m_PVerts[j];
		   pEnd   = pHEdge->pVertEnd;

		   pairP_PHE.first = *pStart;                      //�԰�ߵ������ΪKeyֵ
		   pairP_PHE.second = pHEdge;

		   if (mmapPHE.empty() != true)                    
		   {
			   itVStartH = mmapPHE.lower_bound(*pEnd);       //ָ��Keyֵ��pHEdge�յ���Ϊ���ĵ�һ��pair
			   if (itVStartH != mmapPHE.end())
			   {
			   	   itupperH = mmapPHE.upper_bound(*pEnd);     //ָ��Keyֵ��pHEdge�յ���Ϊ�������һ��pair�ĺ�һ��
				   while(itVStartH != itupperH)
				   {
					   if(itVStartH->second->pHEdgePair ==NULL)//����ð�߻�û�ҵ������
					   {
						   if(itVStartH->second->pVertEnd == pStart)  //�ҵ��յ����pHEdge���İ��Ϊ������
						   {
							   pHEdge->pHEdgePair = itVStartH->second;
							   itVStartH->second->pHEdgePair = pHEdge;  //����߻�Ϊ�����
							   mmapPHE.erase(itVStartH);                  //ɾ�����ҵ�����ߵ�map,2013.6.25��
							   break;
						   }
					   }
					   itVStartH++;
				   }
			   }
		   }
		   m_vecPHEdge.push_back(pHE);             //*****��߱�  
		   if (pHEdge->pHEdgePair == NULL)   //����û�ҵ�����ߣ��ͰѰ����ӵ�map�Եȴ��һ���ߣ�2013.6.25��            
		   {
			   mmapPHE.insert(pairP_PHE);
		   //*****�����߱�
			   PEDGE pEdge = NULL;
			   pEdge = new EDGE();

			   PVERT pPreHEndV,pHEndV;
			   pPreHEndV = pTri->m_PVerts[j];         //���i ��ǰ����ߵ��յ��� ���i �����
			   pHEndV    = pHE->pVertEnd;
			   if (*pPreHEndV < *pHEndV)                                 //�Ƚϰ�ߵ������յ㣬
			   {   pEdge->pVertStart = pPreHEndV;                //*����С�ĵ���Ϊ�ߵ���㣬
				   pEdge->pVertEnd   = pHEndV;                   //*�ϴ����Ϊ�ߵ��յ�
				} 
				else
				{	pEdge->pVertStart = pHEndV;
				    pEdge->pVertEnd   = pPreHEndV;	
				 }
				 pEdge->pHEdgeAdj  = pHE;            //*����ǰ�İ����Ϊ�����ڵ�һ�����

				   //pairP_PE.first = *pEdge->pVertStart; //�Աߵ������ΪKeyֵ
				   //pairP_PE.second = pEdge;

				   //   if (mmapPE.empty() != true)
				   //   {
				   //   itVStart = mmapPE.lower_bound(pairP_PE.first);     //ָ��KeyֵΪ���ĵ�һ��pair
				   //   if (itVStart != mmapPE.end())
				   //   {
				   //      itupper = mmapPE.upper_bound(pairP_PE.first);   //ָ��KeyֵΪ�������һ��pair�ĺ�һ��
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
				   m_vecPEdge.push_back(pEdge);       //*****�߱�	
		   }
		 }
		 pTri->pHEdge = pHE;                                            //*��������Ƭ�����һ����߸���Ƭ	
		 m_vecPFacetTri.push_back(pTri);   //*****���			
		} 
	}

	fclose(file);    //2013.8.15,��Ҫ���ļ��رգ�һ��Ҫע�⣬���������ļ���ռ���ڴ�й¶�����´η����ļ�ʱ��������
	file = NULL;     //2013.8.15,��Ҫ���ļ�ָ��ָ��գ����������ֹ��������ָ�룬��������������ɲ���Ҫ���鷳
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

	//long t2 = GetTickCount();         //�������к�ϵͳʱ��
	//strTime.Format("time:%dms",t2-t1);//ǰ��ʱ���Ϊ��������ʱ��
	//AfxMessageBox(strTime);           //����Ϣ����ʾ
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