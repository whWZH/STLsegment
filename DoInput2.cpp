// DoInput2.cpp : 实现文件
//

#include "stdafx.h"
#include "WZHViewer.h"
#include "DoInput2.h"
#include "afxdialogex.h"


// DoInput2 对话框

IMPLEMENT_DYNAMIC(DoInput2, CDialogEx)

DoInput2::DoInput2(CWnd* pParent /*=NULL*/)
	: CDialogEx(DoInput2::IDD, pParent)
	, theANG(60)
{

}

DoInput2::~DoInput2()
{
}

void DoInput2::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_EDIT1, theANG);
}


BEGIN_MESSAGE_MAP(DoInput2, CDialogEx)
END_MESSAGE_MAP()


// DoInput2 消息处理程序
