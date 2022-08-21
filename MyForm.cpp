#include "MyForm.h"

using namespace System;
using namespace System::Windows::Forms;

[STAThread]
void main(array<String^>^ arg) {
	Application::EnableVisualStyles();
	Application::SetCompatibleTextRenderingDefault(false);
	Project::MyForm form;
	for (int i = 0; i < 5; i++) {
		nowplaying[i] = false;
	}
	Application::Run(%form);
	return;
}