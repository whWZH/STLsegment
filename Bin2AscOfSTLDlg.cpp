// Bin2AscOfSTLDlg.cpp : implementation file
//

#include "stdafx.h"
#include "Bin2AscOfSTL.h"
#include "Bin2AscOfSTLDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

struct FXYZ{float x,y,z;};
struct FACET{FXYZ n;FXYZ dot[3];WORD b;};		//法向量，三节点，保留
/////////////////////////////////////////////////////////////////////////////
// CBin2AscOfSTLDlg dialog

CBin2AscOfSTLDlg::CBin2AscOfSTLDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CBin2AscOfSTLDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(CBin2AscOfSTLDlg)
		// NOTE: the ClassWizard will add member initialization here
	//}}AFX_DATA_INIT
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CBin2AscOfSTLDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CBin2AscOfSTLDlg)
	DDX_Control(pDX, IDC_LIST1, m_lst);
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CBin2AscOfSTLDlg, CDialog)
	//{{AFX_MSG_MAP(CBin2AscOfSTLDlg)
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	//ON_BN_CLICKED(IDC_BUTTON1, OnButton1)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CBin2AscOfSTLDlg message handlers

BOOL CBin2AscOfSTLDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon
	
	// TODO: Add extra initialization here
	
	return TRUE;  // return TRUE  unless you set the focus to a control
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CBin2AscOfSTLDlg::OnPaint() 
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, (WPARAM) dc.GetSafeHdc(), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
}

HCURSOR CBin2AscOfSTLDlg::OnQueryDragIcon()
{
	return (HCURSOR) m_hIcon;
}

void CBin2AscOfSTLDlg::OnButton1() 
{
	// TODO: Add your control notification handler code here
	char st0[10]="solid";
	char st1[20]="facet normal";
	char st2[20]="outer loop";
	char st3[16]="endloop";
	char st4[16]="endfacet";
	char st5[16]="vertex";
	char lf[2];
	char st[80];
	FILE *fp,*fq;
	FACET fct;
	CString DIR,DIR1,DIR2,str;
	int i;
	long l,l1,fcnt,k;
	CFileDialog aDlg(true);
	aDlg.DoModal();
	DIR=aDlg.GetPathName();
	l=DIR.GetLength();
	l1=(aDlg.GetFileName()).GetLength();
	sprintf(st,"%s l=%d l1=%d",DIR,l,l1);//m_lst.AddString(st);
	DIR1=DIR.Left(l-l1);//m_lst.AddString(DIR1);
	DIR2=DIR.Right(l1);//m_lst.AddString(DIR2);
	sprintf(st,"%s",DIR);//m_lst.AddString(st);

	lf[0]=13;lf[1]=10;
//	sprintf(st,"%d",sizeof(FACET));m_lst.AddString(st);return;


	fp = fopen(st,"rb");

	if (!fp){m_lst.AddString("文件打开出错！");return;}
	{
		fseek(fp,0l,SEEK_END); 
		l= ftell(fp);
		l1=(l-84)/sizeof(FACET);
		fseek(fp,80l,SEEK_SET);
		fread(&fcnt,1,sizeof(long),fp);
//		sprintf(st,"l1=%d fcnt=%d",l1,fcnt);m_lst.AddString(st);
		if(l1!=fcnt){m_lst.AddString("二进制ＳＴＬ文件错！");return;}

		sprintf(st,"%sA_%s",DIR1,DIR2);m_lst.AddString(st);
		fq = fopen(st,"wb");sprintf(st,"STL Bianary Format");
		l=strlen(st0);//st0[l]=13;st0[l+1]=10;
		fwrite(st0,1,l,fq);fwrite(lf,1,2,fq);	
		m_lst.AddString(st0);

		for(k=0;k<fcnt;k++)
		{
			fread(&fct,1,sizeof(FACET),fp);
			//normal
			str=st1;
			sprintf(st," %f %f %f",fct.n.x,fct.n.y,fct.n.z);
			str+=st;sprintf(st,"%s",str);
			l=strlen(st);//st[l]=13;st[l+1]=10;
			fwrite(st,1,l,fq);fwrite(lf,1,2,fq);
			m_lst.AddString(st);
			//outer loop
			l=strlen(st2);//st2[l]=13;st2[l+1]=10;
			fwrite(st2,1,l,fq);fwrite(lf,1,2,fq);
			m_lst.AddString(st2);
			//vertex
			for(i=0;i<3;i++)
			{
				str=st5;
				sprintf(st," %f %f %f",fct.dot[i].x,fct.dot[i].y,fct.dot[i].z);
				str+=st;sprintf(st,"%s",str);
				l=strlen(st);//st[l]=13;st[l+1]=10;
				fwrite(st,1,l,fq);fwrite(lf,1,2,fq);
				m_lst.AddString(st);
			}
			//endloop
			l=strlen(st3);//st3[l]=13;st3[l+1]=10;
			fwrite(st3,1,l,fq);fwrite(lf,1,2,fq);
			m_lst.AddString(st3);
			//endfacet
			l=strlen(st4);//st4[l]=13;st4[l+1]=10;
			fwrite(st4,1,l,fq);fwrite(lf,1,2,fq);
			m_lst.AddString(st4);
		}
		fclose(fp);
		fclose(fq);
	}	
}
