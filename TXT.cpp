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
		double DE,K1,K2,SI,Cd,KH,KG;//分别代表是否为边界点，最大曲率K1，最小曲率K2，形状指数SI,曲度Cd
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
	CString s("生成完毕!");
	AfxMessageBox(s);
	system("E:\workplace\ml.py");
}
void CTXT::inputTXT(CSTLModel* pSTLModel)
{
	//////////////////////////////////////////////////////以下为MFC打开文件的对话框
	// CString FilePathName;  CFileDialog dlg(TRUE);///TRUE为OPEN对话框，FALSE为SAVE AS对话框  
	// if(dlg.DoModal()==IDOK)  
	//FilePathName=dlg.GetPathName();
	// ///////////////////////////////////////////////////////
	// /////////////////////////////////////////////////step 1
	// //1. 建立数据类型
	// CvMLData cvml;
	// //2. 读取数据
	//cvml.read_csv(FilePathName);
	// ////3. 指出哪一列代表边界
	// cvml.set_response_idx(0);
	// //4. 选出前4000行数据
	// CvTrainTestSplit cvtts(30, true);
	// //5. Assign the division to the data
	// cvml.set_train_test_split(&cvtts);
	// //6. 声明分类器
	// CvBoost boost;
	// //7. 设置分类器参数，依次为：BOOST分类器类型，分类器个数，样本权值阈值，深度，是否使用代理分裂，
	// boost.train(&cvml, CvBoostParams(CvBoost::REAL, 100, 0, 1, false, 0), false);
	// // 8. 声明一对容器来保存预测结果
	//// vector<float> train_responses, test_responses;
	//// // 9.计算训练误差
	////float fl1 = boost.calc_error(&cvml,CV_TRAIN_ERROR,&train_responses);
	//// // 10. 计算测试误差
 ////   float fl2 = boost.calc_error(&cvml,CV_TEST_ERROR,&test_responses);
 //  boost.save("e:\\trained_boost.xml", "boost");
 //  CvMLData pcvml;
 //  Mat A,B;
 //  pcvml.read_csv("e:\\predict.csv");
 //  pcvml.set_response_idx(0);
 //  A=pcvml.get_values();
 //  B=boost.predict(A);
 //  B=pcvml.get_test_sample_idx();
	/////////////////////读取CSV
	 CString FilePathName;  CFileDialog dlg(TRUE);///TRUE为OPEN对话框，FALSE为SAVE AS对话框  
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