program Project3;

uses
  Vcl.Forms,
  LargeIntegers in '..\..\Arithmetic\LargeIntegers.pas',
  VLargeIntegers in '..\..\Arithmetic\VLargeIntegers.pas',
  GeneralTypes in '..\..\Tools\GeneralTypes.pas',
  HashFunctions in '..\..\Tools\HashFunctions.pas',
  TreeView1 in '..\..\Tools\TreeView1.pas',
  Fp2Arithmetic in '..\..\Extension Fields\Fp2Arithmetic.pas',
  Fp3Arithmetic in '..\..\Extension Fields\Fp3Arithmetic.pas',
  FP4Arithmetic in '..\..\Extension Fields\FP4Arithmetic.pas',
  Fp6Arithmetic in '..\..\Extension Fields\Fp6Arithmetic.pas',
  FP8Arithmetic in '..\..\Extension Fields\FP8Arithmetic.pas',
  FP9Arithmetic in '..\..\Extension Fields\FP9Arithmetic.pas',
  Fp12Arithmetic in '..\..\Extension Fields\Fp12Arithmetic.pas',
  Fp16Arithmetic in '..\..\Extension Fields\Fp16Arithmetic.pas',
  Fp18Arithmetic in '..\..\Extension Fields\Fp18Arithmetic.pas',
  Fp24Arithmetic in '..\..\Extension Fields\Fp24Arithmetic.pas',
  Fp27Arithmetic in '..\..\Extension Fields\Fp27Arithmetic.pas',
  Fp36Arithmetic in '..\..\Extension Fields\Fp36Arithmetic.pas',
  BLS12Curves in '..\..\Curves\BLS12Curves.pas',
  BLS24Curves in '..\..\Curves\BLS24Curves.pas',
  BLS27Curves in '..\..\Curves\BLS27Curves.pas',
  BNCurves in '..\..\Curves\BNCurves.pas',
  ECCFp in '..\..\Curves\ECCFp.pas',
  ECCFp2 in '..\..\Curves\ECCFp2.pas',
  ECCFp3 in '..\..\Curves\ECCFp3.pas',
  ECCFp4 in '..\..\Curves\ECCFp4.pas',
  ECCFp6 in '..\..\Curves\ECCFp6.pas',
  ECCFp9 in '..\..\Curves\ECCFp9.pas',
  KSS16Curves in '..\..\Curves\KSS16Curves.pas',
  KSS18Curves in '..\..\Curves\KSS18Curves.pas',
  KSS36Curves in '..\..\Curves\KSS36Curves.pas',
  MNTCurves in '..\..\Curves\MNTCurves.pas',
  SuperSingularCurves in '..\..\Curves\SuperSingularCurves.pas',
  Unit2 in 'Unit2.pas' {Form2},
  Unit1 in 'Unit1.pas' {Form1},
  Unit3 in 'Unit3.pas' {Form3},
  Unit4 in 'Unit4.pas' {Form4};

{$R *.res}

begin
  Application.Initialize;
  Application.MainFormOnTaskbar := True;
  Application.CreateForm(TForm2, Form2);
  Application.CreateForm(TForm1, Form1);
  Application.CreateForm(TForm2, Form2);
  Application.CreateForm(TForm3, Form3);
  Application.CreateForm(TForm1, Form1);
  Application.CreateForm(TForm3, Form3);
  Application.CreateForm(TForm1, Form1);
  Application.CreateForm(TForm3, Form3);
  Application.CreateForm(TForm4, Form4);
  Application.Run;
end.
