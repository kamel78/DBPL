unit MNTCurves;

interface

uses Vcl.ComCtrls, System.SysUtils, LargeIntegers, Fp12Arithmetic,
  Fp6Arithmetic, Fp3Arithmetic, GeneralTypes, System.classes, Vcl.forms, ECCFp,
  ECCFp3, TreeView1;

Type

  MNTCurvesPairingAlgos = (mnTate, mnAte);
  StandardMNTCurves = (scCurveA, scCurveB, scCurveC);
  GTMNT = Fp6Int;
  G1MNT = FpPoint;
  G2MNT = Fp3Point;
  LargeInt=Lint;
    // https://eprint.iacr.org/2004/165.pdf

  MNTCurvesParamsDefinition = record
    SecurityLevel: String;
    Identifier: String;
    u: String; // the paramater of generation for the MNT curve
    Beta: Integer; // non-square elements of the irredictible polynomial on FP2
    Sigma: string;
    // non-square non-cube elements of the irredictible polynomial on FP3
    A, B: String; // parametres of the curve   y^2=x^3+Ax+B
    TwistMode: TTwistModel;
    GenratorX: String;
    TwistGeneratorBasePointX: String;
    TwistGeneratorSeed: Word;
  end;

  ListOfBNParams = array of MNTCurvesParamsDefinition;

  TMNTCurve = class(TCurve)
    { Definition of a Curve with respect the the Weistrass model   Y^2=X^3+A*X+B (mod P) }
  private
    _A, _B: Lint;
    _u: Lint;
    _LoopMode: LoopPoweringMode; // Loopmode :Negative/Binary Representation
    _FP6PoweringMode: FpPoweringMode;
    // Final Exponentiation Powering Mode : Beuchat Approach, Karabina Approach, Classical Square and Multiply
    _PairingAlgorithm: MNTCurvesPairingAlgos;
    // Pairing Algorithme OptAte,R-Ate.....
    Loop: PLIntArrayForm;
    params:StandardMNTCurves;
    procedure SetA(Value: String);
    function getA: String;
    procedure SetB(Value: String);
    function getB: String;
    procedure Setu(Value: String);
    function getu: String;
    function getRtw: String;
    function getAtw: String;
    function getBtw: String;
    procedure SetSigma(Value: String);
    function getSigma: String;
    procedure SetBeta(Value: Integer);
    function getBeta: Integer;
    function getR: string;
    function getN: string;
    function getTr: string;
    function getLp: string;
    function getHtw: string;
    function getH: string;
    function getP: string;
    procedure setCoord(Value: CoordinatesSystem);
    function getCoord: CoordinatesSystem;
    function GetTwm: TTwistModel;
    procedure RecomputeParametres;
    Procedure FinalPowerOptimalAteMNT(ff: Fp6Int; u: Lint;
      PoweringMode: FpPoweringMode; var Result: Fp6Int);
    // Compute f^((p^18-1)/r)
    Function tatePairing(Pt: FpPoint; Qt: Fp3Point): Fp6Int;
    // Compute Optimal Ate Pairings
    Function AtePairing(Pt: FpPoint; Qt: Fp3Point): Fp6Int;
    procedure SetLoopMode(Value: LoopPoweringMode);
    procedure SetPArams(const Value: StandardMNTCurves);
  public
    // CurveParams:PtrCurveParams;
    Identifier: String;
    /// A symbolic identifier of the curve
    property CoordinatesSystem: CoordinatesSystem read getCoord write setCoord;
    property Beta: Integer read getBeta write SetBeta;
    property Sigma: string read getSigma write SetSigma;
    property u: string read getu write Setu;
    property A: string read getA write SetA;
    property B: string read getB write SetB;
    property Atw: string read getAtw;
    property Btw: string read getBtw;
    property P: string read getP;
    property H: string read getH;
    property N: string read getN;
    property R: String read getR;
    property Rtw: string read getRtw;
    property Htw: string read getHtw;
    property Tr: string read getTr;
    property Lp: string read getLp;
    property TwistMode: TTwistModel read GetTwm;
    constructor Create(AOwner : TComponent); override;
    destructor Destroy;
    procedure SetStandardCurveParametres(inParametres: StandardMNTCurves);
    procedure SetCustomCurveParametres(inParametres: MNTCurvesParamsDefinition);
    procedure GenerateParamsTreeView(Tree: TTreeView); override;
    function GetRandomG1Point: FpPoint;
    function GetRandomG2Point: Fp3Point;
    function GetDefautG1Generator: FpPoint;
    function GetDefautG2Generator: Fp3Point;
    function HashToG1Point(id: string): FpPoint;
    function HashToG2Point(id: String): Fp3Point;
    function Paire(P: FpPoint; Q: Fp3Point): Fp6Int;
    published
    property Fp6PoweringMode: FpPoweringMode read _FP6PoweringMode
      write _FP6PoweringMode;
    property PairingAlgorithm: MNTCurvesPairingAlgos read _PairingAlgorithm
      write _PairingAlgorithm;
    property LoopMode: LoopPoweringMode read _LoopMode write SetLoopMode;
    property Parametres:StandardMNTCurves read params write Setparams;

  end;

const
  // MNT_A_Params_80:MNTCurvesParamsDefinition=(SecurityLevel:'80bits';Identifier:'1';u:'-0x6942ED067F7817837C09';Beta:-2;sigma:'0+u*1';A:'-3';B:'0x77479D33943B5B1F590B54258B72F316B3261D45';GenratorX:'-6');
  MNT_A_Params_80: MNTCurvesParamsDefinition = (SecurityLevel: '80bits';
    Identifier: '1'; u: '-0x6942ED067F7817837C09'; Beta: - 2; Sigma: '0+u*1';
    A: '-3'; B: '0x77479D33943B5B1F590B54258B72F316B3261D45'; GenratorX: '-6');
  MNT_B_Params_106: MNTCurvesParamsDefinition = (SecurityLevel: '106bit';
    Identifier: '1'; u: '-0x7D8EAFA1079E4D49A5A1593976ABC6DF'; Beta: - 6;
    Sigma: '0+u*1'; A: '-3';
    B: '0x6E974D68EF44F266AE3DD5D1F97C497C1D5452D1B074A6C06A25D4E5819CCD1C';
    GenratorX: '-1');
  MNT_C_Params_115: MNTCurvesParamsDefinition = (SecurityLevel: '115bit';
    Identifier: '1'; u: '0x138DD8E76B168FA64C1D73FB49F23E2FA7D9F2B'; Beta: - 5;
    Sigma: '0+u*1'; A: '-3';
    B: '0x05607CD7395B5F49C34A289E4072C37A56601B69C8F64F6BA3F827C87DEE8279BC2E640F16C279';
    GenratorX: '-5');

  MNTParamsList:array[0..2] of string=('MNT_A_Params_80','MNT_B_Params_106','MNT_C_Params_115');
  MNTImplementedPairingAlgos:array[0..1] of string=('Optimal Ate Pairing','Tate');
    procedure Register;

procedure GenerateMNTParametres(Params: MNTCurvesParamsDefinition;
  var Result: PtrCurveParams);
function MNTCurvesParamsDefinitionToSrtingList
  (Params: MNTCurvesParamsDefinition): TstringList;

implementation

{ Generate parameters of a BN curves with precomputed constants from a compact definition }

procedure GenerateMNTParametres(Params: MNTCurvesParamsDefinition;
  var Result: PtrCurveParams);
var
  localP: Lint;
  tmp, tmp1: Lint;
  i: Integer;
  e: Lint;
  twtmp, twtmp1: Fp3Int;
  rho: Lint; // Size of The Miller Loop for Eta Pairings
begin
  if Result = nil then
  begin
    new(Result);
    new(Result.LoopNaf);
    new(Result.LoopBin);
    new(Result.LoopRateBin);
    new(Result.LoopRateNaf);
  end;
  with Result^ do
  begin
    SecurityLevel := Params.SecurityLevel;
    Family := cfMNT;
    if (Pos('0x', Params.u) <> 0) or (Pos('$', Params.u) <> 0) then
      u := _hex_To_LInt(Params.u)
    else
      u := _str_To_LInt(Params.u);
    if Params.Identifier = '1' then
      localP := 4 * u.Sqr + 1
    else
      localP := u.Sqr + 1;

    if not IsLIntPrime(localP) then
      raise Exception.Create
        ('Paramétre t invalide pour la construction de la courbe MNT...')
    else
    begin
      P := localP;
      new(P.InverseFordivision);
      P.Limit := 2 * P.Data.i16[-2];
      P.InverseFordivision^ := (Lint(1) shl (P.Limit * 32)) / P + 1;
      P.id := 'Y';
      R := localP - 2 * u;
      if not IsLIntPrime(R) then
        R := R + 1;

      Tr := 2 * u + 1;
      N := P + 1 - Tr;
      Lp := N;
    end;
    A := _str_To_LInt(Params.A);
    if Params.B[2] = 'x' then
      B := _hex_To_LInt(Params.B)
    else
      B := _str_To_LInt(Params.B);
    Htw := 16 * u.Sqr * u.Sqr + 8 * u.Sqr * u + 12 * u.Sqr;
    Rtw := R;
    H := N / R;
    FieldParam := InitFieldParamsFp3(Params.Beta, P);
    new(FieldParam.MontgomeryData);
    InitMontgomeryStruct(P, FieldParam.MontgomeryData);
    e := Params.Beta;
    if e.IsASqrMod(P, FieldParam.MontgomeryData) then
      Exception.Create(inttostr(Params.Beta) +
        ' is a quadratic Residue modulo P......');
    new(TowerParam);
    TowerParam.pmod8 := P mod 8;
    TowerParam.Sigma.SetFormString(Params.Sigma);
    TowerParam^.FieldParam := FieldParam;
    TowerParam^.Sigma.SetFieldParams(TowerParam^.FieldParam);
    BtwFp3.SetFieldParams(TowerParam^.FieldParam);
    AtwFp3.SetFieldParams(TowerParam^.FieldParam);
    BtwFp3.SetToZero;
    BtwFp3.A := (Params.Beta * Params.Beta * Params.Beta * B) mod P;
    AtwFp3.SetToZero;
    AtwFp3.A := (Params.Beta * Params.Beta * A) mod P;
    LoopBin^ := Lp.ToIntArray;
    LoopNaf^ := Lp.ToNafArray;
    LoopRateBin^ := (2 * u).Absolute.ToIntArray;
    LoopRateNaf^ := (2 * u).Absolute.ToNafArray;
    TwistMode := twMType;
    /// ///
    if (Pos('0x', Params.GenratorX) <> 0) or (Pos('$', Params.GenratorX) <> 0)
    then
      BasePointX := _hex_To_LInt(Params.GenratorX)
    else
      BasePointX := _str_To_LInt(Params.GenratorX);
    /// Constructing the generator for G1
    _Sqr_LInt(BasePointX, tmp);
    _Add_LInt(tmp, A, tmp);
    _Mul_LInt(tmp, BasePointX, tmp1);
    _Add_LInt(tmp1, B, tmp);
    _Mod_LInt(tmp, P, tmp);
    if tmp.IsASqrMod(P, FieldParam.MontgomeryData) then
    begin
      ModSquareLInt(tmp, P, BasePointY, FieldParam.MontgomeryData, false);
      BasePointY := P - BasePointY;
    end
    else
      raise Exception.Create
        ('Invalide parameter for the generator of the curve');
    /// Constructing the generator for G2
    TwistBasePointX_Fp3.SetFieldParams(TowerParam^.FieldParam);
    TwistBasePointy_Fp3.SetFieldParams(TowerParam^.FieldParam);

    TwistBasePointX_Fp3.A.SetToRandom(TowerParam^.FieldParam.P,
      Params.TwistGeneratorSeed);
    TwistBasePointX_Fp3.B.SetToRandom(TowerParam^.FieldParam.P,
      Params.TwistGeneratorSeed);
    TwistBasePointX_Fp3.c.SetToRandom(TowerParam^.FieldParam.P,
      Params.TwistGeneratorSeed);

    repeat
      _Inc_LInt(TwistBasePointX_Fp3.A, 1);
      _Sqr_FP3(TwistBasePointX_Fp3, twtmp);
      _Add_FP3(twtmp, AtwFp3, twtmp);
      _Mul_FP3(twtmp, TwistBasePointX_Fp3, twtmp1);
      _Add_FP3(twtmp1, BtwFp3, twtmp);
    until twtmp.IsASquare;
    _Sqrt_FP3(twtmp, TwistBasePointy_Fp3);
    _Neg_FP3(TwistBasePointy_Fp3, twtmp1);
    if _Compare_FP3(TwistBasePointy_Fp3, twtmp1) = 1 then
      TwistBasePointy_Fp3 := twtmp1;
    TowerParam^.FrobeniusP_Const[0].SetFieldParams(TowerParam^.FieldParam);
    TowerParam^.FrobeniusP_Const[0] := TowerParam^.Sigma.Pow((P - 1) / 3);
    TowerParam^.FrobeniusP_Const[1].SetFieldParams(TowerParam^.FieldParam);
    TowerParam^.FrobeniusP_Const[1] := TowerParam^.FrobeniusP_Const[0] *
      TowerParam^.FrobeniusP_Const[0];
    u := 2 * u;
  end;
end;

function MNTCurvesParamsDefinitionToSrtingList
  (Params: MNTCurvesParamsDefinition): TstringList;
begin
  Result := TstringList.Create;
  Result.Add('Parameter u :' + Params.u);
  if Params.Beta < 0 then
    Result.Add('Irrudictible Polynomial  for Fp2: X^2-' +
      inttostr(Abs(Params.Beta)))
  else
    Result.Add('Irrudictible Polynomial  for Fp2: X^2+' +
      inttostr(Abs(Params.Beta)));
  Result.Add('Irrudictible Polynomial  for Fp6: X^3+' + Params.Sigma);
  if Params.B[1] <> '-' then
    Result.Add('Curve : y^2=x^3+' + Params.B)
  else
    Result.Add('Curve : y^2=x^3' + Params.B);
  if Params.TwistMode = twDType then
    Result.Add('Twiste Mode: D-Type')
  else
    Result.Add('Twiste Mode: M-Type');
  Result.Add('Hamming Weight of u (Binary)' +
    inttostr(HammingWeight(Lint(Params.u), false)));
  Result.Add('Hamming Weight of u (NAF)' +
    inttostr(HammingWeight(Lint(Params.u), True)));
end;

{ TMNTCurve }

procedure TMNTCurve.GenerateParamsTreeView(Tree: TTreeView);
var
  s: string;
  R: real;
  G: FpPoint;
  Gtw: Fp3Point;
  tmp: TTreeNode;
begin
  Tree.Items.Clear;
  tmp := Tree.Items.Add(nil, 'Curve');
  Tree.Items.AddChild(tmp, 'Family : MNT');
  if CurveParams.A = 1 then
    s := 'y^2=x^3+x'
  else
    s := 'y^2=' + CurveParams.A.ToDecimalString + '*x^3+x';
  if not(_IsNull(CurveParams.B)) then
    s := s + CurveParams.B.ToDecimalString;
  Tree.Items.AddChild(tmp, 'Equation : ' + s);
  tmp := Tree.Items.Add(nil, 'Initial Parameter x0');
  Tree.Items.AddChild(tmp, 'x0= ' + CurveParams.u.ToHexString);
  Tree.Items.AddChild(tmp, 'Hamming Weight (Bin) =' +
    inttostr(HammingWeight(CurveParams.u, false)));
  Tree.Items.AddChild(tmp, 'Hamming Weight (Neg) =' +
    inttostr(HammingWeight(CurveParams.u, True)));
  tmp := Tree.Items.Add(nil, 'Field Fp');
  Tree.Items.AddChild(tmp, 'Prime P = ' + CurveParams.P.ToHexString);
  Tree.Items.AddChild(tmp, 'Size of Fp : ' + inttostr(CurveParams.P.BitLength)
    + ' bit');
  tmp := Tree.Items.Add(nil, 'Field E(Fp)');
  Tree.Items.AddChild(tmp, 'Order (#E(Fp)) N = ' + CurveParams.N.ToHexString);
  Tree.Items.AddChild(tmp, 'Frobenius Trace Tr : ' +
    CurveParams.Tr.ToHexString);
  tmp := Tree.Items.Add(nil, 'Base Fields G1 : E(Fp)[R]');
  Tree.Items.AddChild(tmp, 'Order of G1 R = ' + CurveParams.R.ToHexString);
  Tree.Items.AddChild(tmp, 'Size of G1 : ' + inttostr(CurveParams.R.BitLength)
    + ' bit');
  Tree.Items.AddChild(tmp, 'Cofactor of G1 H = ' + CurveParams.H.ToHexString);
  tmp := Tree.Items.Add(nil, 'Twist');
  Tree.Items.AddChild(tmp, 'Equation of the Twist : ' + 'y^2=x^3+' +
    CurveParams.AtwFp3.ToHexString + '*x+' + CurveParams.BtwFp3.ToHexString);
  Tree.Items.AddChild(tmp, 'Degree of the Twist : Quadratic');
  s := 'M-Type';
  Tree.Items.AddChild(tmp, 'Type of the Twist : ' + s);
  tmp := Tree.Items.Add(nil, 'Field G2 (Twist) :E''(Fp3)[R]');
  Tree.Items.AddChild(tmp, 'Order of G2 R = ' + CurveParams.Rtw.ToHexString);
  Tree.Items.AddChild(tmp, 'Size of G2 : ' + inttostr(CurveParams.Rtw.BitLength)
    + ' bit');
  Tree.Items.AddChild(tmp, 'Cofactor of G2 Htw = ' +
    CurveParams.Htw.ToHexString);
  tmp := Tree.Items.Add(nil, 'Targted Field (GT):');
  Tree.Items.AddChild(tmp, 'Extension Field : Fp6');
  Tree.Items.AddChild(tmp, 'Size of GT : ' +
    inttostr(6 * CurveParams.P.BitLength) + ' bit');
  tmp := Tree.Items.Add(nil, 'Security Level');
  Tree.Items.AddChild(tmp, 'Security in G1 : ' +
    inttostr(CurveParams.R.BitLength div 2) + ' bit').ImageIndex := 2;
  Tree.Items.AddChild(tmp, 'Security in G2 : ' +
    inttostr(CurveParams.R.BitLength div 2) + ' bit').ImageIndex := 2;
  R := 1.92 * exp((1 / 3) * ln(6 * CurveParams.P.BitLength)) *
    sqr(exp((1 / 3) * (ln(ln(6 * CurveParams.P.BitLength) / ln(2)))));
  Tree.Items.AddChild(tmp, 'Security in GT (GNFS) : ' + inttostr(Round(R)) +
    ' bit').ImageIndex := 2;
  tmp := Tree.Items.Add(nil, 'Implemented Pairings');
  Tree.Items.AddChild(tmp, 'Tate');
  tmp := Tree.Items.Add(nil, 'Tower Construction');
  if CurveParams.FieldParam.Beta < 0 then
    Tree.Items.AddChild(tmp, 'Fp3<u>=ExstensionField<u,|u^3-' +
      inttostr(Abs(CurveParams.FieldParam.Beta)) + '>')
  else
    Tree.Items.AddChild(tmp, 'Fp3<u>=ExstensionField<u,|u^3+' +
      inttostr(CurveParams.FieldParam.Beta) + '>');
  Tree.Items.AddChild(tmp, 'Fp6<v>=ExstensionField<v,|v^2-' +
    CurveParams.TowerParam.Sigma.ToHexString + '>');
  tmp := Tree.Items.Add(nil, ' = ' + floattostrf(CurveParams.P.BitLength /
    CurveParams.R.BitLength, ffGeneral, 2, 4));
  tmp := Tree.Items.Add(nil, 'Default G1 Generator');
  G.SetCurveParams(CurveParams);
  G.SetToDefaultGenerator;
  Gtw.SetCurveParams(CurveParams);
  Gtw.SetToDefaultGenerator;
  Tree.Items.AddChild(tmp, G.ToHexString);
  tmp := Tree.Items.Add(nil, 'Default G2 Generator');
  Tree.Items.AddChild(tmp, Gtw.ToHexString);
  Tree.FullExpand;
  Tree.Selected := Tree.Items[0];
  //Tree.SetFocus;
end;


{ ******************************************************************************* }
constructor TMNTCurve.Create(AOwner : TComponent);
begin
  inherited create(Aowner);
  CurveParams := nil;
  SetStandardCurveParametres(scCurveA);
  PairingAlgorithm := mnAte;
  CoordinatesSystem := csAffine;
  Fp6PoweringMode := pmNormal;
  SetLoopMode(_LoopMode);
end;

{ ******************************************************************************* }
destructor TMNTCurve.Destroy;
begin
  Dispose(CurveParams);
  inherited;
end;

{ ******************************************************************************* }
function TMNTCurve.getA: String;
begin
  Result := _A.ToHexString;
end;

{ ******************************************************************************* }
function TMNTCurve.getAtw: String;
begin
  Result := CurveParams^.AtwFp3.ToHexString;
end;

{ ******************************************************************************* }
function TMNTCurve.getB: String;
begin
  Result := _B.ToHexString;
end;

{ ******************************************************************************* }
function TMNTCurve.getBeta: Integer;
begin
  Result := CurveParams^.FieldParam.Beta;
end;

{ ******************************************************************************* }
function TMNTCurve.getBtw: String;
begin
  Result := CurveParams^.BtwFp3.ToHexString;
end;

{ ******************************************************************************* }
function TMNTCurve.getCoord: CoordinatesSystem;
begin
  Result := CurveParams^.CoordSys;
end;

{ ******************************************************************************* }
function TMNTCurve.GetDefautG1Generator: FpPoint;
begin
  Result.SetToDefaultGenerator;
end;

{ ******************************************************************************* }
function TMNTCurve.GetDefautG2Generator: Fp3Point;
begin
  Result.SetToDefaultGenerator;
end;

{ ******************************************************************************* }
function TMNTCurve.getH: string;
begin
  Result := CurveParams^.H.ToHexString;
end;

{ ******************************************************************************* }
function TMNTCurve.getHtw: string;
begin
  Result := CurveParams^.Htw.ToHexString;
end;

{ ******************************************************************************* }
function TMNTCurve.getLp: string;
begin
  Result := CurveParams^.Lp.ToHexString;
end;

{ ******************************************************************************* }
function TMNTCurve.getN: string;
begin
  Result := CurveParams^.N.ToHexString;
end;

{ ******************************************************************************* }
function TMNTCurve.getP: string;
begin
  Result := CurveParams^.P.ToHexString;
end;

{ ******************************************************************************* }
function TMNTCurve.getR: string;
begin
  Result := CurveParams^.R.ToHexString;
end;

function TMNTCurve.GetRandomG1Point: FpPoint;
begin
  Result.SetCurveParams(CurveParams);
  Result.SetAsRandomTorsionPoint;
end;

function TMNTCurve.GetRandomG2Point: Fp3Point;
begin
  Result.SetCurveParams(CurveParams);
  Result.SetAsRandomTorsionPoint;
end;

{ ******************************************************************************* }
function TMNTCurve.getRtw: String;
begin
  Result := CurveParams^.Rtw.ToHexString;
end;

{ ******************************************************************************* }
function TMNTCurve.getSigma: String;
begin
  Result := CurveParams^.TowerParam.Sigma.ToHexString;
end;

{ ******************************************************************************* }
function TMNTCurve.getTr: string;
begin
  Result := CurveParams^.Tr.ToHexString;
end;

{ ******************************************************************************* }
function TMNTCurve.GetTwm: TTwistModel;
begin
  Result := CurveParams.TwistMode;
end;

{ ******************************************************************************* }
function TMNTCurve.getu: String;
begin
  Result := _u.ToHexString;
end;

{ ******************************************************************************* }
function TMNTCurve.HashToG1Point(id: string): FpPoint;
begin
  Result.SetCurveParams(CurveParams);
  Result.SetAsTorsionFromString(id);
end;

{ ******************************************************************************* }
function TMNTCurve.HashToG2Point(id: String): Fp3Point;
begin
  Result.SetCurveParams(CurveParams);
  Result.SetAsTorsionFromString(id);
end;

{ ******************************************************************************* }
procedure TMNTCurve.RecomputeParametres;
var
  tmp: MNTCurvesParamsDefinition;
begin
  tmp.u := _u.ToHexString;
  tmp.Beta := Beta;
  tmp.Sigma := Sigma;
  tmp.A := _A.ToHexString;
  tmp.B := _B.ToHexString;
  tmp.Identifier := Identifier;
  GenerateMNTParametres(tmp, CurveParams);
end;

{ ******************************************************************************* }
procedure TMNTCurve.SetStandardCurveParametres(inParametres: StandardMNTCurves);
begin
  case inParametres of
    scCurveA:
      begin
        GenerateMNTParametres(MNT_A_Params_80, CurveParams);
        _A := _str_To_LInt(MNT_A_Params_80.A);
        _B := _hex_To_LInt(MNT_A_Params_80.B);
        _u := _hex_To_LInt(MNT_A_Params_80.u);
      end;
    scCurveB:
      begin
        GenerateMNTParametres(MNT_B_Params_106, CurveParams);
        _A := _str_To_LInt(MNT_B_Params_106.A);
        _B := _hex_To_LInt(MNT_B_Params_106.B);
        _u := _hex_To_LInt(MNT_B_Params_106.u);
      end;
    scCurveC:
      begin
        GenerateMNTParametres(MNT_C_Params_115, CurveParams);
        _A := _str_To_LInt(MNT_C_Params_115.A);
        _B := _hex_To_LInt(MNT_C_Params_115.B);
        _u := _hex_To_LInt(MNT_C_Params_115.u);
      end;
  end;
end;

{ ******************************************************************************* }
procedure TMNTCurve.SetA(Value: String);
begin
  RecomputeParametres;
end;

{ ******************************************************************************* }
procedure TMNTCurve.SetB(Value: String);
begin
  RecomputeParametres;
end;

{ ******************************************************************************* }
procedure TMNTCurve.SetBeta(Value: Integer);
begin
  RecomputeParametres;
end;

{ ******************************************************************************* }
procedure TMNTCurve.setCoord(Value: CoordinatesSystem);
begin
  CurveParams^.CoordSys := Value;
end;

{ ******************************************************************************* }
procedure TMNTCurve.SetCustomCurveParametres(inParametres
  : MNTCurvesParamsDefinition);
begin
  GenerateMNTParametres(inParametres, CurveParams);
  _A := _str_To_LInt(inParametres.A);
  _B := _str_To_LInt(inParametres.B);
  _u := _hex_To_LInt(inParametres.u);
end;

{ ******************************************************************************* }
procedure TMNTCurve.SetSigma(Value: String);
begin
  RecomputeParametres;
end;

{ ******************************************************************************* }
procedure TMNTCurve.Setu(Value: String);
begin
  RecomputeParametres;
end;

{ ********* Pairing function :Compute the pairing according to the specified algorithm*********** }
function TMNTCurve.Paire(P: FpPoint; Q: Fp3Point): Fp6Int;
begin
  if (P.Infinity) or (Q.Infinity) then
    raise Exception.Create('Can''t compute pairing for infinit points ....');
  case PairingAlgorithm of
    mnTate:
      Result := tatePairing(P, Q);
    mnAte:
      Result := AtePairing(P, Q);
  end;
end;

{ ******************** Set the loop mode of the pairing :Binary / Negatve Representation************* }
procedure TMNTCurve.SetLoopMode(Value: LoopPoweringMode);
begin
  _LoopMode := Value;
  case Value of
    lpmBinary:
      begin
        case PairingAlgorithm of
          mnTate:
            Loop := CurveParams.LoopBin;
          mnAte:
            Loop := CurveParams.LoopRateBin;
        end;

      end;
    lpmNaf:
      begin
        case PairingAlgorithm of
          mnTate:
            Loop := CurveParams.LoopNaf;
          mnAte:
            Loop := CurveParams.LoopRateNaf;
        end;
      end;
  end;
end;

procedure TMNTCurve.SetPArams(const Value: StandardMNTCurves);
var pr:TComponent;
    i:integer;
begin
SetStandardCurveParametres(Value);
pr:=Self.Owner;
for i:=0 to pr.ComponentCount-1 do begin
                                   if (pr.Components[i] is TTreeView1)and (TTreeView1(pr.Components[i]).Curve = self)
                                   then GenerateParamsTreeView(TTreeView(pr.Components[i]));
                                   end;
params:=value;
end;

{ ******************** Final Exponentiation Step :Compute f^((p^6-1)/r) *********************** }
Procedure TMNTCurve.FinalPowerOptimalAteMNT(ff: Fp6Int; u: Lint;
  PoweringMode: FpPoweringMode; var Result: Fp6Int); // Compute f^((p^18-1)/r)
var
  tmp, tmp1, invf: Fp6Int;
begin
  _Pow_FP6_P(ff, tmp);
  _Mul_FP6(ff, tmp, tmp);
  _Pow_FP6_P_i(tmp, 3, tmp1);
  _Inv_FP6(tmp, tmp);
  _Mul_FP6(tmp1, tmp, tmp);
  _Pow_FP6_P(tmp, tmp1);
  _Pow_FP6(tmp, u.Absolute, tmp);
  if _IsNeg(u) then
    _Inv_FP6(tmp, tmp);
  _Mul_FP6(tmp, tmp1, Result);
end;

{ ******************** Compute the Tate Pairing *********************** }
Function TMNTCurve.tatePairing(Pt: FpPoint; Qt: Fp3Point): Fp6Int;
var
  T, mPt, tmpPt: FpPoint;
  i: Integer;
  f: Fp6Int;
  u: Fp3Int;
begin
  f.SetTowerParams(CurveParams.TowerParam);
  f.SetToOne;
  T := Pt;
  T.ComputeLineAtQPFp6 := True;
  T.LineAtQPFp6.SetTowerParams(CurveParams.TowerParam);
  T.SetMNTPairingPointCoordinates(Qt.X, Qt.Y);
  // Untwist Qt
  u.SetFieldParams(CurveParams.FieldParam);
  u.SetToZero;
  u.B := 1;
  T.CurrentXQMNT := Lint(CurveParams.FieldParam.Beta)
    .InversModulo(CurveParams.P) * T.CurrentXQMNT;
  T.CurrentYQMNT := Lint(CurveParams.FieldParam.Beta)
    .InversModulo(CurveParams.P).Sqr * T.CurrentYQMNT * u;
  //
  _Neg_Fp_Point(Pt, mPt);
  for i := Length(Loop^) - 2 Downto 0 do
  begin
    case CoordinatesSystem of
      csAffine:
        _Double_Affine_Fp_Point(T, T, false);
      csProjective:
        _Double_projective_Fp_Point(T, T, false);
    end;
    f := f * f;
    f := f * T.LineAtQPFp6;
    if (Loop^[i] <> 0) then
    begin
      if Loop^[i] = 1 then
        tmpPt := Pt
      else
        tmpPt := mPt;
      case CoordinatesSystem of
        csAffine:
          _Add_Affine_Fp_Point(tmpPt, T, T, false);
        csProjective:
          _Add_projective_Fp_Point(T, tmpPt, T, false);
      end;
      f := f * T.LineAtQPFp6;
    end;
  end;
  FinalPowerOptimalAteMNT(f, CurveParams.u, Fp6PoweringMode, Result);
end;

{ ******************** Compute the Ate Pairing *********************** }
Function TMNTCurve.AtePairing(Pt: FpPoint; Qt: Fp3Point): Fp6Int;
var
  T, mQt, tmpQt: Fp3Point;
  i: Integer;
  f: Fp6Int;
  u: Fp3Int;
begin
  f.SetTowerParams(CurveParams.TowerParam);
  f.SetToOne;
  T := Qt;
  T.ComputeLigneAtFp6 := True;
  T.ComputeLigneValue := false;
  T.LineAtQPFp6.SetTowerParams(CurveParams.TowerParam);
  T.SetPairingPointCoordinates(Pt.X * Abs(CurveParams.FieldParam.Beta),
    Pt.Y * Abs(CurveParams.FieldParam.Beta));
  _Neg_Fp3_Point(Qt, mQt);
  for i := Length(Loop^) - 2 Downto 0 do
  begin
    case CoordinatesSystem of
      csAffine:
        _Double_Affine_Fp3_Point(T, T, false);
      csProjective:
        _Double_Affine_Fp3_Point(T, T, false);
      /// Projective ccordinates for MNT Curves has not been implemented (For the Ate pairing)
      /// since the same pointes are used for KSS18, and has not the same projective formulation
      /// (the MNT curves has non nule A parameters) So implemented Projective curves for MNT in the same
      /// unit will make the code very  harder to understand, In addition the gain with projective coordinates for MNT is not important (1 ms)
    end;
    f := f * f;
    f := f * T.LineAtQPFp6;
    if (Loop^[i] <> 0) then
    begin
      if Loop^[i] = 1 then
        tmpQt := Qt
      else
        tmpQt := mQt;
      case CoordinatesSystem of
        csAffine:
          _Add_Affine_Fp3_Point(tmpQt, T, T, false);
        csProjective:
          _Add_Affine_Fp3_Point(tmpQt, T, T, false);
        /// Projective ccordinates for MNT Curves has not been implemented (For the Ate pairing)
        /// since the same pointes are used for KSS18, and has not the same projective formulation
        /// (the MNT curves has non nule A parameters) So implemented Projective curves for MNT in the same
        /// unit will make the code very  harder to understand, In addition the gain with projective coordinates for MNT is not important (1 ms)
      end;
      f := f * T.LineAtQPFp6;
    end;
  end;
  FinalPowerOptimalAteMNT(f, CurveParams.u, Fp6PoweringMode, Result);
end;

procedure Register;
begin
  RegisterComponents('Pairings', [TMNTCurve]);
end;

end.
