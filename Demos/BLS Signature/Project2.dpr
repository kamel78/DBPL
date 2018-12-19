program Project2;

uses
  Vcl.Forms,
  Unit2 in 'Unit2.pas' {Form2},
  LargeIntegers in '..\..\My Github pairings\Arithmetic\LargeIntegers.pas',
  VLargeIntegers in '..\..\My Github pairings\Arithmetic\VLargeIntegers.pas',
  BLS12Curves in '..\..\My Github pairings\Curves\BLS12Curves.pas',
  BLS24Curves in '..\..\My Github pairings\Curves\BLS24Curves.pas',
  BLS27Curves in '..\..\My Github pairings\Curves\BLS27Curves.pas',
  BNCurves in '..\..\My Github pairings\Curves\BNCurves.pas',
  ECCFp in '..\..\My Github pairings\Curves\ECCFp.pas',
  ECCFp2 in '..\..\My Github pairings\Curves\ECCFp2.pas',
  ECCFp3 in '..\..\My Github pairings\Curves\ECCFp3.pas',
  ECCFp4 in '..\..\My Github pairings\Curves\ECCFp4.pas',
  ECCFp6 in '..\..\My Github pairings\Curves\ECCFp6.pas',
  ECCFp9 in '..\..\My Github pairings\Curves\ECCFp9.pas',
  KSS16Curves in '..\..\My Github pairings\Curves\KSS16Curves.pas',
  KSS18Curves in '..\..\My Github pairings\Curves\KSS18Curves.pas',
  KSS36Curves in '..\..\My Github pairings\Curves\KSS36Curves.pas',
  MNTCurves in '..\..\My Github pairings\Curves\MNTCurves.pas',
  SuperSingularCurves in '..\..\My Github pairings\Curves\SuperSingularCurves.pas',
  Fp2Arithmetic in '..\..\My Github pairings\Extension Fields\Fp2Arithmetic.pas',
  Fp3Arithmetic in '..\..\My Github pairings\Extension Fields\Fp3Arithmetic.pas',
  FP4Arithmetic in '..\..\My Github pairings\Extension Fields\FP4Arithmetic.pas',
  Fp6Arithmetic in '..\..\My Github pairings\Extension Fields\Fp6Arithmetic.pas',
  FP8Arithmetic in '..\..\My Github pairings\Extension Fields\FP8Arithmetic.pas',
  FP9Arithmetic in '..\..\My Github pairings\Extension Fields\FP9Arithmetic.pas',
  Fp12Arithmetic in '..\..\My Github pairings\Extension Fields\Fp12Arithmetic.pas',
  Fp16Arithmetic in '..\..\My Github pairings\Extension Fields\Fp16Arithmetic.pas',
  Fp18Arithmetic in '..\..\My Github pairings\Extension Fields\Fp18Arithmetic.pas',
  Fp24Arithmetic in '..\..\My Github pairings\Extension Fields\Fp24Arithmetic.pas',
  Fp27Arithmetic in '..\..\My Github pairings\Extension Fields\Fp27Arithmetic.pas',
  Fp36Arithmetic in '..\..\My Github pairings\Extension Fields\Fp36Arithmetic.pas',
  GeneralTypes in '..\..\My Github pairings\Tools\GeneralTypes.pas',
  HashFunctions in '..\..\My Github pairings\Tools\HashFunctions.pas',
  TreeView1 in '..\..\My Github pairings\Tools\TreeView1.pas';

{$R *.res}

begin
  Application.Initialize;
  Application.MainFormOnTaskbar := True;
  Application.CreateForm(TForm2, Form2);
  Application.Run;
end.
