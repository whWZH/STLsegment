// Bin2AscOfSTL.cpp : Defines the class behaviors for the application.
//

#include "stdafx.h"
#include "Bin2AscOfSTL.h"
#include "Bin2AscOfSTLDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CBin2AscOfSTLApp

BEGIN_MESSAGE_MAP(CBin2AscOfSTLApp, CWinApp)
	//{{AFX_MSG_MAP(CBin2AscOfSTLApp)
	//}}AFX_MSG
	ON_COMMAND(ID_HELP, CWinApp::OnHelp)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CBin2AscOfSTLApp construction

CBin2AscOfSTLApp::CBin2AscOfSTLApp()
{
}

/////////////////////////////////////////////////////////////////////////////
// The one and only CBin2AscOfSTLApp object

CBin2AscOfSTLApp theApp;

/////////////////////////////////////////////////////////////////////////////
// CBin2AscOfSTLApp initialization

BOOL CBin2AscOfSTLApp::InitInstance()
{
	// Standard initialization

	CBin2AscOfSTLDlg dlg;
	m_pMainWnd = &dlg;
	int nResponse = dlg.DoModal();
	if (nResponse == IDOK)
	{
	}
	else if (nResponse == IDCANCEL)
	{
	}

	// Since the dialog has been closed, return FALSE so that we exit the
	//  application, rather than start the application's message pump.
	return FALSE;
}
