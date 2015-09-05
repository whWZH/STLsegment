#include "stdafx.h"
#include "TXT.h"
#include "WZHViewerDoc.h"
#include "GLView.h"
#include "GeomBase.h"
#include "solid.h"
#include <math.h>
CTXT::CTXT()
{	

}
CTXT::~CTXT()
{

}
void CTXT::outputTXT(CSTLModel* pSTLModel)
{
	ofstream ocout;
	/*vector<double> sdf_gram;
	vector<vector<PVERT>> RBLOOP;
	pCSDF->SDF_gram(4,RBLOOP,sdf_gram,pSTLModel->m_vecPFacetTri,pSTLModel->vec_box,pSTLModel->m_vecPVert);*/
	ocout.open("e:\\theResult.csv");
	for (int i=0;i<pSTLModel->m_vecPVert.size();i++)
	{
		double DE,K1,K2,SI,Cd,KH,KG;//�ֱ�����Ƿ�Ϊ�߽�㣬�������K1����С����K2����״ָ��SI,����Cd
		PVERT P=pSTLModel->m_vecPVert[i];
		KG=Ccalcubase::CalcuVerGausCurvature(P);
		KH = Ccalcubase::CalcuVerMeanCurvature(P);
		Ccalcubase::CalcuVerPrinCurvature(P, K1, K2);
		SI = 2 / (3.1415926*atan(K1 / K2));
		Cd = sqrt(0.5*(K1*K1 + K2*K2));
		float DE1,K11,K21,SI1,Cd1,KH1,KG1;
		DE=0;
		DE1=DE;K11=K1;K21=K2;SI1=SI;Cd1=Cd;KH1=KH;KG1=KG;
		if (P->bused==1)
		{
			DE1=1;
			ocout<<1<<","<<K11<<","<<K21<<","<<SI1<<","<<Cd1<<","<<KH1<<","<<KG1<<"\n";
		}
		else
		{
			DE1=0;
			ocout<<0<<","<<K11<<","<<K21<<","<<SI1<<","<<Cd1<<","<<KH1<<","<<KG1<<"\n";
		}
		//ocout<<DE1<<","<<K11<<","<<K21<<","<<SI1<<","<<Cd1<<","<<KH1<<","<<KG1<<"\n";
	}
	/*for (int i=0;i<sdf_gram.size();i++)
	{
		int num=sdf_gram[i];
		ocout<<num<<"\n";
	}*/
	ocout.close();
	CString s("�������!");
	AfxMessageBox(s);
	system("E:\workplace\ml.py");
}
void CTXT::inputTXT(CSTLModel* pSTLModel)
{
	//////////////////////////////////////////////////////����ΪMFC���ļ��ĶԻ���
	// CString FilePathName;  CFileDialog dlg(TRUE);///TRUEΪOPEN�Ի���FALSEΪSAVE AS�Ի���  
	// if(dlg.DoModal()==IDOK)  
	//FilePathName=dlg.GetPathName();
	// ///////////////////////////////////////////////////////
	// /////////////////////////////////////////////////step 1
	// //1. ������������
	// CvMLData cvml;
	// //2. ��ȡ����
	//cvml.read_csv(FilePathName);
	// ////3. ָ����һ�д���߽�
	// cvml.set_response_idx(0);
	// //4. ѡ��ǰ4000������
	// CvTrainTestSplit cvtts(30, true);
	// //5. Assign the division to the data
	// cvml.set_train_test_split(&cvtts);
	// //6. ����������
	// CvBoost boost;
	// //7. ���÷���������������Ϊ��BOOST���������ͣ�����������������Ȩֵ��ֵ����ȣ��Ƿ�ʹ�ô�����ѣ�
	// boost.train(&cvml, CvBoostParams(CvBoost::REAL, 100, 0, 1, false, 0), false);
	// // 8. ����һ������������Ԥ����
	//// vector<float> train_responses, test_responses;
	//// // 9.����ѵ�����
	////float fl1 = boost.calc_error(&cvml,CV_TRAIN_ERROR,&train_responses);
	//// // 10. ����������
 ////   float fl2 = boost.calc_error(&cvml,CV_TEST_ERROR,&test_responses);
 //  boost.save("e:\\trained_boost.xml", "boost");
 //  CvMLData pcvml;
 //  Mat A,B;
 //  pcvml.read_csv("e:\\predict.csv");
 //  pcvml.set_response_idx(0);
 //  A=pcvml.get_values();
 //  B=boost.predict(A);
 //  B=pcvml.get_test_sample_idx();
	/////////////////////��ȡCSV
	 CString FilePathName;  CFileDialog dlg(TRUE);///TRUEΪOPEN�Ի���FALSEΪSAVE AS�Ի���  
	 if(dlg.DoModal()==IDOK)  
	FilePathName=dlg.GetPathName();
	 ifstream fstrm(FilePathName);
	 string line;
	 vector < int >  IDS;
	 while (getline(fstrm, line))
	 {
		 int ID;
		 sscanf(line.c_str(), "%d", &ID);
		 IDS.push_back(ID);
	 };
	 for (int i = 0; i < pSTLModel->m_vecPVert.size();i++)
	 {
		 if (IDS[i]==1)
		 {
			 pSTLModel->m_vecPVert[i]->bused = 1;
		 }
	 }
}