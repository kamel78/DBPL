unit BLS12Curves;

interface

uses  Vcl.ComCtrls,VCL.dialogs, System.SysUtils,VLargeIntegers,LargeIntegers,Fp12Arithmetic,Fp6Arithmetic,Fp4Arithmetic,Fp2Arithmetic,
      GeneralTypes,System.classes,ECCFp2,ECCFp, Treeview1;


Type

   BLS12CurvesPairingAlgos=(bls12pOptAte);
   G1BLS12=FpPoint;
   G2BLS12=Fp2Point;
   GTBLS12=Fp12Int;
   LargeInt=Lint;
   StandardBLS12Curves=(sc128_1,sc128_2_raz,sc128_3,sc128_4);
   BLS12CurvesParamsDefinition=record
                            SecurityLevel:String;
                            u:String;  // the paramater of generation for the BLS12 curve
                            Beta:Integer; // non-square elements of the irredictible polynomial on FP2
                            Sigma:string; // non-square non-cube elements of the irredictible polynomial on FP4
                            Gamma:string; // non-square non-cube elements of the irredictible polynomial on FP8
                            A,B:String; // parametres of the curve   y^2=x^3+Ax+B
                            TwistMode:TTwistModel;
                            GenratorX:String;
                            TwistGeneratorSeed:Word;
                            end;
   ListOfBLS12Params=array of BLS12CurvesParamsDefinition;

 TBLS12CurvePairing=class (TCurve)
                      {Definition of a Curve with respect the the Weistrass model   Y^2=X^3+A*X+B (mod P)}
              private
              _A,_B:Lint;
              _u:Lint;
              _LoopMode:LoopPoweringMode;                     //  Loopmode :Negative/Binary Representation
              _Fp12PoweringMode:FpPoweringMode;             // Final Exponentiation Powering Mode : Beuchat Approach, Karabina Approach, Classical Square and Multiply
              _PairingAlgorithm:BLS12CurvesPairingAlgos;         //  Pairing Algorithme OptAte,R-Ate.....
              Loop:PLIntArrayForm;
              params:StandardBLS12Curves;
              procedure SetA(Value:String);
              function getA:String;
              procedure SetB(Value:String);
              function getB:String;
              procedure Setu(Value:String);
              function getu:String;
              function getRtw:String;
              function getAtw:String;
              function getBtw:String;
              procedure SetSigma(Value:String);
              function getSigma:String;
              function getGamma:String;
              procedure SetGamma(Value:String);
              procedure SetBeta(Value:integer);
              function getBeta:integer;
              function getR:string;
              function getN:string;
              function getTr:string;
              function getLp:string;
              function getHtw:string;
              function getH:string;
              function getP:string;

              procedure setCoord(value:CoordinatesSystem);
              function getCoord:CoordinatesSystem;
              function GetTwm:TTwistModel;
              procedure RecomputeParametres;
              procedure SetPArams(const Value: StandardBLS12Curves);

              Procedure FinalPowerOptimalAteBLS12(ff:Fp12Int;u:LInt;PoweringMode:FpPoweringMode;var Result:Fp12Int); // Compute f^((p^12-1)/r)
              Function OptAtePairing(Pt:FpPoint;Qt:Fp2Point):Fp12Int;  // Compute Optimal Ate Pairings
              procedure SetLoopMode(Value:LoopPoweringMode);

            public
              SecurityLevel:String;    /// A symbolic identifier of the curve
              PerferedPoweringMode:LoopPoweringMode;

              property Beta:integer read getBeta write SetBeta;
              property Sigma:string read Getsigma write SetSigma;
              property Gamma:string read GetGamma Write SetGamma;
              property u :string read Getu write Setu;
              property A:string read getA write setA;
              property B:string read getB write setB;
              property Atw :string read GetAtw;
              property Btw :string read GetBtw;
              property P:string read GetP;
              property H:string read GetH;
              property N:string read getN;
              property R:String read getR;
              property Rtw:string read getRtw;
              property Htw:string read getHtw;
              property Tr:string read getTr;
              property Lp:string read getLp;
              property TwistMode:TTwistModel read GetTwm;

              constructor Create(AOwner : TComponent); override;
              destructor Destroy;
              procedure SetStandardCurveParametres(inParametres:StandardBLS12Curves);
              procedure SetCustomCurveParametres(inParametres:BLS12CurvesParamsDefinition);
              procedure GenerateParamsTreeView(Tree:TTreeView);override;
              function GetRandomG1Point:FpPoint;
              function GetRandomG2Point:Fp2Point;
              function GetDefautG1Generator:FpPoint;
              function GetDefautG2Generator:Fp2Point;
              function HashToG1Point(id:string):FpPoint;
              function HashToG2Point(id:String):Fp2Point;

              function Paire(P:FpPoint;Q:Fp2Point):Fp12Int;
              published
              property Fp12PoweringMode:FpPoweringMode read _Fp12PoweringMode write _Fp12PoweringMode;
              property PairingAlgorithm:BLS12CurvesPairingAlgos read _PairingAlgorithm write _PairingAlgorithm ;
              property LoopMode:LoopPoweringMode read _LoopMode write SetLoopMode;
              property CoordinatesSystem:CoordinatesSystem read getCoord write SetCoord;
              property Parametres:StandardBLS12Curves read params write Setparams;

            end;

const

  BLS12_128_1_Params:BLS12CurvesParamsDefinition=(SecurityLevel:'128bit';u:'0xC000000000040405';Beta:-2;sigma:'0+u*1';A:'0';B:'1';TwistMode:twDType;GenratorX:'6';TwistGeneratorSeed:55);
  BLS12_128_3_Params:BLS12CurvesParamsDefinition=(SecurityLevel:'128bit';u:'-0xC00000000003D799';Beta:-1;sigma:'1+u*1';A:'0';B:'1';TwistMode:twDType;GenratorX:'4';TwistGeneratorSeed:55);
  BLS12_Razvan_128_Params:BLS12CurvesParamsDefinition=(SecurityLevel:'128bit';u:'-0x1FFFFFFBFFFE00000000';Beta:-1;sigma:'1+u*1';A:'0';B:'4';TwistMode:twMType;GenratorX:'5';TwistGeneratorSeed:55);
  BLS12_128_4_Params:BLS12CurvesParamsDefinition=(SecurityLevel:'128bit';u:'0x4000000000000438';Beta:-1;sigma:'1+u*2';A:'0';B:'6';TwistMode:twDType;GenratorX:'1';TwistGeneratorSeed:55);

  BLS12ParamsList:array[0..3] of string=('BLS12_128_1_Params','BLS12_Razvan_128_Params','BLS12_128_3_Params','BLS12_128_4_Params');
  BLS12ImplementedPairingAlgos:array[0..0] of string=('Optimal Ate Pairing');

  procedure ComputeBLS12Parametres(Params: BLS12CurvesParamsDefinition; var Result:PtrCurveParams);

  procedure Register;

implementation


Procedure FastLoadParams(param:StandardBLS12Curves;var Result:PtrCurveParams);
var rho,tmp,e:Lint;
    i:integer;
begin
if Result=nil then begin
                   new(Result);
                   new(Result.LoopBin);
                   new(Result.LoopNaf);
                   new(Result.LoopRateBin);
                   new(Result.LoopRateNaf);
                   new(Result.LoopEtaNaf);
                   new(Result.LoopEtaBin);
                   new(Result.LoopTateNaf);
                   new(Result.LoopTateBin);
                   end;
with Result^ do begin
case param of
sc128_1:begin
               SecurityLevel:='128bit';
               Family:=cfBLS12;
               u:='0xC000000000040405';
               P:='0xF3000000001E7E83700001982F49B3247B62097A5541F34E632E14502A27D57F2C73B1B4FD2335009F8C7CCB2009B35';
               new(p.InverseFordivision);
               p.Limit:=2*p.Data.i16[-2];
               p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
               p.id:='Y';
               r:='0x510000000006C6C8700000366CBD8753D0C24516665180738B131967B0FFAA59';
               Tr:='0xC000000000040406';
               Lp:= '0xC000000000040405';
               n:='0xF3000000001E7E83700001982F49B3247B62097A5541F34E632E14502A27D57F2C73B1B4FD2334F49F8C7CCB1FC9730';
               rho:='0x9366C48000000005B696800000000013A700000000000017';
               htw:='0x2D90000000079FA08B00008EDD4DEFA2CAF9DDD4F05E28F1982477EACA1000BFED27E9E9F945D20525B0F5DC317201B552111E134E7F55D4C7B53D5A515C691';
               rtw:=r;
               A:=0;
               B:=1;
               H:='0x30000000000202020000000560100AB0';
               FieldParam:=InitFieldParams(-2,p);
               New(FieldParam.MontgomeryData);
               InitMontgomeryStruct(P,FieldParam.MontgomeryData);
               New(TowerParam);
               TowerParam.Sigma.SetFormString('0+u*1');
               Atw.SetFormString('0');
               TowerParam^.FieldParam:=FieldParam;
               Atw.SetFieldParams(TowerParam^.FieldParam);
               Atw.SetFormString('0');
               TowerParam^.Sigma.SetFieldParams(TowerParam^.FieldParam);
               Btw.SetFieldParams(TowerParam^.FieldParam);
               TwistMode:=twDType;
               Btw:=B*TowerParam.Sigma.Inverse;
               LoopBin^:=Lp.ToIntArray;
               LoopNaf^:=Lp.ToNafArray;
               LoopRateBin^:=(tr-1).ToIntArray;
               LoopRateNaf^:=(tr-1).ToNafArray;
               LoopEtaBin^:=rho.ToIntArray;
               LoopEtaNaf^:=rho.ToNafArray;
               LoopTateBin^:=n.ToIntArray;
               LoopTateNaf^:=n.ToNafArray;
               BasePointX:='6';
               BasePointY:='0x6A3633EAAA5B86F3B6B0F457B21573F3D57FD5B07C79FB2612F4335AEE1E4D0257FF76D54CDF3EFE91F72A3DC5D7AB5';
               TwistBasePointX.SetFieldParams(TowerParam^.FieldParam);
               TwistBasePointy.SetFieldParams(TowerParam^.FieldParam);
               TwistBasePointX.SetFormString('0xA810FD57AAD2C7D51BC7EA21EBB9A0CECE0F370D96BC1054B6B58E14BEF0872391D13AB2F5F6DC40F08117F8E8F825+u*0xA810FD57AAD2C7D51BC7EA21EBB9A0CECE0F370D96BC1054B6B58E14BEF0872391D13AB2F5F6DC40F08117F8E8F81F');
               TwistBasePointY.SetFormString('0x1EC1D7426303DDB5E470D77434FF2519BB85BF28056F556E455BA279257FDC950E6C44E50110F764AF31E8BC719EDDA+u*0xED6F0EC1957C3848815BA954F845867B830ED5A08B3280F2927CF306CE16BB398D1BCA03523C5348CE2A2C22F87C600');
        end;
sc128_2_raz:begin
               SecurityLevel:='128bit';
               Family:=cfBLS12;
               u:='-0x1FFFFFFBFFFE00000000';
               P:='0x15555545554D5A555A55D69414935FBD6F1E32D8BACCA47B14848B42A8DFFA5C1CC00F26AA91557F00400020000555554AAAAAAC0000AAAAAAAB';
               new(p.InverseFordivision);
               p.Limit:=2*p.Data.i16[-2];
               p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
               p.id:='Y';
               r:='0xFFFFFF7FFFC0180017FE05FD000E801FC017FFC80001100007FEFFFEFFFFC0000000000000001';
               Tr:='-0x1FFFFFFBFFFDFFFFFFFF';
               Lp:= '0x1FFFFFFBFFFE00000000';
               n:='0x15555545554D5A555A55D69414935FBD6F1E32D8BACCA47B14848B42A8DFFA5C1CC00F26AA91557F00400020000555556AAAAAA7FFFEAAAAAAAB';
               rho:='0x8E521A9886E502411DC1AF70120000017E80600000000001';
               htw:='0x1C71C6FFFFF1D38E4555CA343641384EF449EF40F2574F227721F5F081128BA285EBA2329EBA89241E66E72C7A130E799C48DC91839D531BC31A9B720071E755538E31D5538E371C70E38E38E5';
               rtw:=r;
               A:=0;
               B:=4;
               H:='0x1555554FFFFD55AAAB01556AAA7FFFEAAAAAAAB';
               FieldParam:=InitFieldParams(-1,p);
               New(FieldParam.MontgomeryData);
               InitMontgomeryStruct(P,FieldParam.MontgomeryData);
               New(TowerParam);
               TowerParam.Sigma.SetFormString('1+u*1');
               Atw.SetFormString('0');
               TowerParam^.FieldParam:=FieldParam;
               Atw.SetFieldParams(TowerParam^.FieldParam);
               Atw.SetFormString('0');
               TowerParam^.Sigma.SetFieldParams(TowerParam^.FieldParam);
               Btw.SetFieldParams(TowerParam^.FieldParam);
               TwistMode:=twMType;
               Btw:=B*TowerParam.Sigma;
               LoopBin^:=Lp.ToIntArray;
               LoopNaf^:=Lp.ToNafArray;
               LoopRateBin^:=(tr-1).ToIntArray;
               LoopRateNaf^:=(tr-1).ToNafArray;
               LoopEtaBin^:=rho.ToIntArray;
               LoopEtaNaf^:=rho.ToNafArray;
               LoopTateBin^:=n.ToIntArray;
               LoopTateNaf^:=n.ToNafArray;
               BasePointX:='5';
               BasePointY:='0x493BA66B9DE5026933364FC307E10BABCA0CE952E1BBB8BA8D3022728905676F4AFAEB0453C46CF80E06101944B8D53BD919602379D77481160';
               TwistBasePointX.SetFieldParams(TowerParam^.FieldParam);
               TwistBasePointy.SetFieldParams(TowerParam^.FieldParam);
               TwistBasePointX.SetFormString('0x982C4C9C050D0FA008CDC243F576689AD5B9178E6F01D66C0A13DEAF7B110C0CA590184DFFEE114EF490B6347D6E1EEE81FBFA270BE0793093D+u*0x982C4C9C050D0FA008CDC243F576689AD5B9178E6F01D66C0A13DEAF7B110C0CA590184DFFEE114EF490B6347D6E1EEE81FBFA270BE0793093A');
               TwistBasePointY.SetFormString('0x745136C3C39150C999B8C77D4368F80AE5AE0D9AFF2ED09D6466E5AEA7790861ED0E2D17C6F8B6AAB9F087A8DB3F0056802DC4C418B5F2C2F87+u*0x1354725EEC21208F09BA08239603DCFA11765F2ADEC72A1E319C7255D039AAEA2F1A75C62C6363C65654C20BC5DB89C711332829DBBE3C0379B7');
          end;
sc128_3:begin
               SecurityLevel:='128bit';
               Family:=cfBLS12;
               u:='-0xC00000000003D799';
               P:='0xF3000000001D2D3460000175AC2E2E3239F8637C111DF63C3EC7A5D698DD5C3F4908DA58BAAEE9117CACCB82892BB33';
               new(p.InverseFordivision);
               p.Limit:=2*p.Data.i16[-2];
               p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
               p.id:='Y';
               r:='0x5100000000067BD230000031D2A82DDCD0AA2897B7A5500F6D3046EB75676A71'  ;
               Tr:='-0xC00000000003D798';
               Lp:='0xC00000000003D799';
               n:='0xF3000000001D2D3460000175AC2E2E3239F8637C111DF63C3EC7A5D698DD5C3F4908DA58BAAEE91D7CACCB8289692CC';
               rho:='0x9003628ECA2A09940C9D88B06C3A2EC9E84970AC7B2A2717';
               htw:='0x2D90000000074B4D69000082C91B8276323C01B5EA79D64664BF96BA5AFD20530E081623687FE4769E33D96CE3399B8D2FA989EA24AB4897C84994B659237C5';
               rtw:=r;
               A:=0;
               B:=1;
               H:='0x300000000001EBCD00000004EBBAAD8C';
               FieldParam:=InitFieldParams(-1,p);
               New(FieldParam.MontgomeryData);
               InitMontgomeryStruct(P,FieldParam.MontgomeryData);
               New(TowerParam);
               TowerParam.Sigma.SetFormString('1+u*1');
               Atw.SetFormString('0');
               TowerParam^.FieldParam:=FieldParam;
               Atw.SetFieldParams(TowerParam^.FieldParam);
               Atw.SetFormString('0');
               TowerParam^.Sigma.SetFieldParams(TowerParam^.FieldParam);
               Btw.SetFieldParams(TowerParam^.FieldParam);
               TwistMode:=twDType;
               Btw:=B*TowerParam.Sigma.Inverse;
               LoopBin^:=Lp.ToIntArray;
               LoopNaf^:=Lp.ToNafArray;
               LoopRateBin^:=(tr-1).ToIntArray;
               LoopRateNaf^:=(tr-1).ToNafArray;
               LoopEtaBin^:=rho.ToIntArray;
               LoopEtaNaf^:=rho.ToNafArray;
               LoopTateBin^:=n.ToIntArray;
               LoopTateNaf^:=n.ToNafArray;
               BasePointX:='4';
               BasePointY:='0x3E46D9E8905FB9E38690AE3A5B06BC82431FBE6CE0FBFF4537FF29B577BD2F76DCBFFCF4CEA9987CCA7FBDA833CF85B';
               TwistBasePointX.SetFieldParams(TowerParam^.FieldParam);
               TwistBasePointy.SetFieldParams(TowerParam^.FieldParam);
               TwistBasePointX.SetFormString('0xA810FD57AAFCF1B71BC7EE724F2A3F16FB43F6D61B3BB29943835D46E83FAF1FFF2C263B44805A254C77410BC4B827+u*0xA810FD57AAFCF1B71BC7EE724F2A3F16FB43F6D61B3BB29943835D46E83FAF1FFF2C263B44805A254C77410BC4B823');
               TwistBasePointY.SetFormString('0x8DF3549AD7E23102CD0AC7B441FEB620A49AE14FFAE39237F70EFF11F182169FBF39E04A8BF4D98F683FE6469A4E993+u*0x1EFB66DD4FBF8E1CA81F088B9170D05F313904567385EBBEC799C0D5DF2CABE333737B4DD074DFBF9CDFC22DDC17AC0');
          end;
sc128_4:begin
               SecurityLevel:='128bit';
               Family:=cfBLS12;
               u:='0x4000000000000438';
               P:='0x5555555555557712AAAAAAAAB0399DAAAAAAAB27B122EAD5555B83321325C5557F079AA30437571FE6C5DCC182013';
               new(p.InverseFordivision);
               p.Limit:=2*p.Data.i16[-2];
               p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
               p.id:='Y';
               r:='0x100000000000043800000000006AC97F00000004B15ABE40000013CC36443C1';
               Tr:='0x4000000000000439';
               Lp:='0x4000000000000438';
               n:='0x5555555555557712AAAAAAAAB0399DAAAAAAAB27B122EAD5555B83321325C5557F079AA30437531FE6C5DCC181BDB';
               rho:='0x90000D86C2F7D9E58CEAC8F763D9687528C67610746ACBEF';
               htw:='0x1C71C71C71C72B700000000003751E08E38E39583119278E38ED2960D956F87248BBAC81DA33CE0C4A3C21AFEE02A993399CEED3ADE8C5297D8C24228205';
               rtw:=r;
               A:=0;
               B:=6;
               H:='0x5555555555556092AAAAAAAAAB0969B';
               FieldParam:=InitFieldParams(-1,p);
               New(FieldParam.MontgomeryData);
               InitMontgomeryStruct(P,FieldParam.MontgomeryData);
               New(TowerParam);
               TowerParam.Sigma.SetFormString('1+u*2');
               Atw.SetFormString('0');
               TowerParam^.FieldParam:=FieldParam;
               Atw.SetFieldParams(TowerParam^.FieldParam);
               Atw.SetFormString('0');
               TowerParam^.Sigma.SetFieldParams(TowerParam^.FieldParam);
               Btw.SetFieldParams(TowerParam^.FieldParam);
               TwistMode:=twDType;
               Btw:=B*TowerParam.Sigma.Inverse;
               LoopBin^:=Lp.ToIntArray;
               LoopNaf^:=Lp.ToNafArray;
               LoopRateBin^:=(tr-1).ToIntArray;
               LoopRateNaf^:=(tr-1).ToNafArray;
               LoopEtaBin^:=rho.ToIntArray;
               LoopEtaNaf^:=rho.ToNafArray;
               LoopTateBin^:=n.ToIntArray;
               LoopTateNaf^:=n.ToNafArray;
               BasePointX:='1';
               BasePointY:='0x16305FDD9C96B63C604EE3CF5EAA3F156060C343EAD2750759132EE4EDD3ACE6F51292B42621DF8CE125C42C53039';
               TwistBasePointX.SetFieldParams(TowerParam^.FieldParam);
               TwistBasePointy.SetFieldParams(TowerParam^.FieldParam);
               TwistBasePointX.SetFormString('0x2BBA8025940FF09F11D727B24FE24838FA58E3A6A38FBF3C4E3B6B5EE002AA4A293CA6FE1FF604C0DC3214497A009+u*0x2BBA8025940FF09F11D727B24FE24838FA58E3A6A38FBF3C4E3B6B5EE002AA4A293CA6FE1FF604C0DC3214497A006');
               TwistBasePointY.SetFormString('0x714365CB775797FE4B359240E6833F796C43B61B3A4E5BB99E7F7D4C56D1117923EF775210DCC94D71C911AB3C4E+u*0x1BBE25D07C993DBADA70E2715A8DD3495A4429679F7BB99C14F4362401C3C5DADAA36984EEBE3AE424B45D279621');
          end;

end;
TowerParam^.FrobeniusP_Const[0].SetFieldParams(TowerParam^.FieldParam);
TowerParam^.FrobeniusP_Const[0]:=TowerParam^.Sigma.Pow((P-1)/6);
for i:=1 to 4 do begin
                 TowerParam^.FrobeniusP_Const[i].SetFieldParams(TowerParam^.FieldParam);
                 TowerParam^.FrobeniusP_Const[i]:=TowerParam^.FrobeniusP_Const[i-1]*TowerParam^.FrobeniusP_Const[0];
                 end;
for i:=0 to 4 do begin
                 TowerParam^.FrobeniusP2_Const[i].SetFieldParams(TowerParam^.FieldParam);
                 TowerParam^.FrobeniusP3_Const[i].SetFieldParams(TowerParam^.FieldParam);
                 TowerParam^.FrobeniusP2_Const[i]:=TowerParam^.FrobeniusP_Const[i]*TowerParam^.FrobeniusP_Const[i].Conjugate;
                 TowerParam^.FrobeniusP3_Const[i]:=TowerParam^.FrobeniusP_Const[i]*TowerParam^.FrobeniusP2_Const[i];
                 end;
for i:=0 to 2 do begin
                 e:=P.Pow(i+1);
                 FrobeniusMapConstX[i].SetFieldParams (TowerParam.FieldParam);
                 case TwistMode of
                                 twDType: begin
                                          FrobeniusMapConstX[i]:=TowerParam^.Sigma.Pow((e-1)/3);
                                          FrobeniusMapConstY[i]:=TowerParam^.Sigma.Pow((e-1)/2);
                                          end;
                                 twMType:begin
                                         FrobeniusMapConstX[i]:=(TowerParam^.Sigma.Pow((e-1)/3)).inverse;
                                         FrobeniusMapConstY[i]:=TowerParam^.Sigma.Pow((e-1)/2).inverse;
                                         end;
                            end;
                 end;
_Pow_LInt(P,power,4);
_Pow_LInt(P,tmp,2);
_Sub_LInt(power,tmp,power);
_Inc_LInt(power,1);
_Div_Mod_LInt(power,R,power,tmp);
end;
end;


procedure ComputeBLS12Parametres(Params: BLS12CurvesParamsDefinition; var Result:PtrCurveParams);
var localP:LInt;
    u6,u5,u3,tmp,tmp1:LInt;
    i:integer;
    m2:HLInt;
    e:Lint;
    twtmp,twtmp1:Fp2Int;
    tau:array[0..4] of HLint;
    s:string;
begin
if Result=nil then begin
                   new(Result);
                   new(Result.LoopBin);
                   new(Result.LoopNaf);
                   new(Result.LoopTateNaf);
                   new(Result.LoopTateBin);
                   end;
with Result^ do begin
                SecurityLevel:=Params.SecurityLevel;
                Family:=cfBLS12;
                if (Pos('0x',Params.u)<>0)or(Pos('$',Params.u)<>0) then u:=_hex_To_LInt(Params.u)
                else u:=_str_To_LInt(Params.u);
                u3:=(u.Sqr)*u;
                u5:=u3*u.Sqr;
                u6:=u5*u;
                LocalP:=(u6-2*u5+2*u3+u+1)/3;
                if not IsLIntPrime(LocalP) then raise Exception.Create('Paramétre t invalide pour la construction de la courbe BLS...')
                else begin
                     p:=Localp;
                     p.Limit:=2*p.Data.i16[-2];
                     new(p.InverseFordivision);
                     p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
                     p.id:='Y';
                     Tr:=(u+1);
                     Lp:=u.Absolute;
                     r:=u*(u3-u)+1;
                     n:=p+1-tr;
                     end;
                A:=_Str_To_LInt(Params.A);
                B:=_Str_To_LInt(Params.B);
                FieldParam:=InitFieldParams(Params.Beta,p);
                e:=params.beta;
                if e.IsASqrMod(P,FieldParam.MontgomeryData) then Exception.Create(inttostr(Params.beta)+' is a quadratic Residue modulo P......');
                New(FieldParam.MontgomeryData);
                InitMontgomeryStruct(P,FieldParam.MontgomeryData);
                //computing number of points on E(FP2)with cofactor
                Tau[0]:=tr.Sqr-2*p;
                Tau[1]:=p.Sqr;
                _Sqrt_HLint(((4*Tau[1]-Tau[0].Sqr)/3),Tau[2]);
                s:=Tau[2].ToHexString;
                m2:=Tau[1]+1-((Tau[0]-3*Tau[2])/2);
                _VHCopy_LInt(m2,tmp1);

                Htw:=tmp1/r;    // Twist Cofactor
                Rtw:=r;         ///  The true order is m2=r*htw , so we can choose a subgourp of order R (the cofactor is htw) (https://eprint.iacr.org/2005/133.pdf)
                H:=n/r;
                New(TowerParam);
                TowerParam^.FieldParam:=FieldParam;
                TowerParam^.Sigma.SetFieldParams(TowerParam^.FieldParam);
                TowerParam.Sigma.SetFormString(Params.Sigma);
                     ///   should test if beta and sigma are valide parametres
                Twtmp.SetFieldParams(FieldParam);
                _Pow_FP2(TowerParam.Sigma,(Result.P.Sqr-1)/2,Twtmp);
                if (Twtmp.a=1)and(Twtmp.b=0) then raise Exception.Create('Paramétre t invalide pour la construction de la courbe BLS...Sigma is a square');
                _Pow_FP2(Result.TowerParam.Sigma,(Result.P.Sqr-1)/3,Twtmp);
                if (Twtmp.a=1)and(Twtmp.b=0) then raise Exception.Create('Paramétre t invalide pour la construction de la courbe BLS...Sigma is a Cube');
                    ///
                TowerParam.pmod8:=p mod 8;
                Btw.SetFieldParams(FieldParam);
                TwistMode:=Params.TwistMode;
                if TwistMode=twMType then _Mul_FP_FP2(B,TowerParam.Sigma,Btw) else //raise Exception.Create('Invalide parameter for the generator of the curve, this code supports only D-Type Twistes ....');
                _Mul_FP_FP2(B,TowerParam.Sigma.Inverse,Btw);
                LoopBin^:=Lp.ToIntArray;
                LoopNaf^:=Lp.ToNafArray;
                LoopTateBin^:=n.ToIntArray;
                LoopTateNaf^:=n.ToNafArray;
                if (Pos('0x',Params.GenratorX)<>0)or(Pos('$',Params.GenratorX)<>0) then BasePointX:=_hex_To_LInt(Params.GenratorX)
                else BasePointX:=_Str_To_LInt(Params.GenratorX);
                ///   Constructing the base point of the generator for G1
                _Sqr_LInt(BasePointX,tmp);
                _Add_LInt(tmp,A,tmp);
                _Mul_LInt(tmp,BasePointX,tmp1);
                _Add_LInt(tmp1,B,tmp);
                _Mod_LInt(tmp,P,tmp);
                if tmp.IsASqrMod(P,FieldParam.MontgomeryData)
                                              then begin
                                                   ModSquareLInt(tmp,P, BasePointY,FieldParam.MontgomeryData,false);
                                                   if P-BasePointY<BasePointY then BasePointY:=P-BasePointY;
                                                   end
                else raise Exception.Create('Invalide parameter for the generator of the curve');
                ///   Constructing the base point for G2 Generator
                TwistBasePointX.SetFieldParams(TowerParam^.FieldParam);
                TwistBasePointy.SetFieldParams(TowerParam^.FieldParam);

                TwistBasePointX.a.SetToRandom(TowerParam^.FieldParam.p,params.TwistGeneratorSeed);
                TwistBasePointX.b.SetToRandom(TowerParam^.FieldParam.p,params.TwistGeneratorSeed);
                repeat
                _Inc_LInt(TwistBasePointX.a,1);
                _Sqr_FP2(TwistBasePointX,twtmp);
                _Add_FP2(twtmp,Atw,twtmp);
                _Mul_FP2(twtmp,TwistBasePointX,twtmp1);
                _Add_FP2(twtmp1,Btw,twtmp);
                until twtmp.IsASquare;
                _Sqrt_FP2(twtmp,TwistBasePointy);
                _Neg_FP2(TwistBasePointy,twtmp1);
                if _Compare_FP2(TwistBasePointy,twtmp1)=1 then TwistBasePointy:=twtmp1;
                TowerParam^.FrobeniusP_Const[0].SetFieldParams(TowerParam^.FieldParam);
                TowerParam^.FrobeniusP_Const[0]:=TowerParam^.Sigma.Pow((P-1)/6);
                for i:=1 to 4 do begin
                                 TowerParam^.FrobeniusP_Const[i].SetFieldParams(TowerParam^.FieldParam);
                                 //if i=3 then showmessage(TowerParam^.FrobeniusP_Const[i-1].toHexString);

                                 _Mul_FP2(TowerParam^.FrobeniusP_Const[i-1],TowerParam^.FrobeniusP_Const[0],TowerParam^.FrobeniusP_Const[i]);
                                 //TowerParam^.FrobeniusP_Const[i]:=TowerParam^.FrobeniusP_Const[i-1]*TowerParam^.FrobeniusP_Const[0];
                                 end;
                             //exit;
                for i:=0 to 4 do begin
                                 TowerParam^.FrobeniusP2_Const[i].SetFieldParams(TowerParam^.FieldParam);
                                 TowerParam^.FrobeniusP3_Const[i].SetFieldParams(TowerParam^.FieldParam);
                                 TowerParam^.FrobeniusP2_Const[i]:=TowerParam^.FrobeniusP_Const[i]*TowerParam^.FrobeniusP_Const[i].Conjugate;
                                 TowerParam^.FrobeniusP3_Const[i]:=TowerParam^.FrobeniusP_Const[i]*TowerParam^.FrobeniusP2_Const[i];
                                 end;
                for i:=0 to 2 do begin
                                 e:=P.Pow(i+1);
                                 FrobeniusMapConstX[i].SetFieldParams (TowerParam.FieldParam);
                                 case TwistMode of
                                 twDType: begin
                                          FrobeniusMapConstX[i]:=TowerParam^.Sigma.Pow((e-1)/3);
                                          FrobeniusMapConstY[i]:=TowerParam^.Sigma.Pow((e-1)/2);
                                          end;
                                 twMType:begin
                                         FrobeniusMapConstX[i]:=(TowerParam^.Sigma.Pow((e-1)/3)).Inverse;///
                                         FrobeniusMapConstY[i]:=(TowerParam^.Sigma.Pow((e-1)/2)).Inverse;   //////
                                         end;
                                 end;
                                 end;
                end;
end;





{ TBLS12Curve }


procedure TBLS12CurvePairing.GenerateParamsTreeView(Tree:TTreeView);
var s:string;
    r:real;
    G:FpPoint;Gtw:Fp2Point;
    tmp:TTreeNode;
begin
Tree.Items.Clear;
tmp:=Tree.Items.Add(nil,'Curve');
Tree.Items.AddChild(tmp,'Family : BLS12');
if CurveParams.B<0 then s:='y^2=x^3-'+CurveParams.B.Absolute.ToDecimalString
else s:='y^2=x^3-'+CurveParams.B.ToDecimalString;
Tree.Items.AddChild(tmp,'Equation : '+s);
tmp:=Tree.Items.Add(nil,'Initial Parameter x0');
Tree.Items.AddChild(tmp,'x0= '+curveparams.u.ToHexString );
Tree.Items.AddChild(tmp,'Hamming Weight (Bin) ='+IntToStr(HammingWeight(CurveParams.u,False)));
Tree.Items.AddChild(tmp,'Hamming Weight (Neg) ='+IntToStr(HammingWeight(CurveParams.u,True)));
tmp:=Tree.Items.Add(nil,'Field Fp');
Tree.Items.AddChild(tmp,'Prime P = '+CurveParams.P.ToHexString );
Tree.Items.AddChild(tmp,'Size of Fp : '+inttostr(CurveParams.P.BitLength)+' bit');
tmp:=Tree.Items.Add(nil,'Field E(Fp)');
Tree.Items.AddChild(tmp,'Order (#E(Fp)) N = '+CurveParams.N.ToHexString );
Tree.Items.AddChild(tmp,'Frobenius Trace Tr : '+CurveParams.Tr.ToHexString);
tmp:=Tree.Items.Add(nil,'Base Fields G1 : E(Fp)[R]');
Tree.Items.AddChild(tmp,'Order of G1 R = '+CurveParams.R.ToHexString);
Tree.Items.AddChild(tmp,'Size of G1 : '+inttostr(CurveParams.R.BitLength)+' bit');
Tree.Items.AddChild(tmp,'Cofactor of G1 H = '+CurveParams.H.ToHexString);
tmp:=Tree.Items.Add(nil,'Twist');
Tree.Items.AddChild(tmp,'Equation of the Twist : '+'y^2=x^3+'+Curveparams.Btw.toHexString);
Tree.Items.AddChild(tmp,'Degree of the Twist : Sextic');
if Curveparams.TwistMode=twMType then S:='M-Type' else S:='D-Type';
Tree.Items.AddChild(tmp,'Type of the Twist : '+S);
tmp:=Tree.Items.Add(nil,'Field G2 (Twist) :E''(Fp2)[R]');
Tree.Items.AddChild(tmp,'Order of G2 R = '+CurveParams.Rtw.ToHexString);
Tree.Items.AddChild(tmp,'Size of G2 : '+inttostr(CurveParams.Rtw.BitLength)+' bit');
Tree.Items.AddChild(tmp,'Cofactor of G2 Htw = '+CurveParams.Htw.ToHexString);
tmp:=Tree.Items.Add(nil,'Targted Field (GT):');
Tree.Items.AddChild(tmp,'Extension Field : Fp12');
Tree.Items.AddChild(tmp,'Size of GT : '+inttostr(12*CurveParams.P.BitLength)+' bit');
tmp:=Tree.Items.Add(nil,'Security Level');
Tree.Items.AddChild(tmp,'Security in G1 : '+inttostr(CurveParams.R.BitLength div 2)+' bit').ImageIndex:=2;
Tree.Items.AddChild(tmp,'Security in G2 : '+inttostr(CurveParams.R.BitLength div 2)+' bit').ImageIndex:=2;
r:=1.92* exp((1/3)*ln(12*CurveParams.P.BitLength))*sqr(exp((1/3)*(ln(ln(12*CurveParams.P.BitLength)/ln(2)))));
Tree.Items.AddChild(tmp,'Security in GT (GNFS) : '+inttostr(Round(r))+' bit').ImageIndex:=2;
tmp:=Tree.Items.Add(nil,'Implemented Pairings');
Tree.Items.AddChild(tmp,'Otp-Ate');
tmp:=Tree.Items.Add(nil,'Tower Construction');
if Curveparams.FieldParam.Beta<0 then Tree.Items.AddChild(tmp,'Fp2<u>=ExstensionField<u,|u^2+'+inttostr(Abs(Curveparams.FieldParam.Beta))+'>')
else Tree.Items.AddChild(tmp,'Fp2<u>=ExstensionField<u,|u^2-'+inttostr(Curveparams.FieldParam.Beta)+'>');
Tree.Items.AddChild(tmp,'Fp6<v>=ExstensionField<v,|v^3-('+Curveparams.TowerParam.Sigma.toHexString+')>');
Tree.Items.AddChild(tmp,'Fp12<w>=ExstensionField<w,|w^2-v>');
tmp:=Tree.Items.Add(nil,' = '+floattostrf(CurveParams.P.BitLength/CurveParams.R.BitLength,ffGeneral, 2, 4 ));
tmp:=Tree.Items.Add(nil,'Default G1 Generator');
G.SetCurveParams (CurveParams);
G.SetToDefaultGenerator;
Gtw.SetCurveParams(CurveParams);
Gtw.SetToDefaultGenerator;
Tree.Items.AddChild(tmp,G.toHexString);
tmp:=Tree.Items.Add(nil,'Default G2 Generator');
Tree.Items.AddChild(tmp,Gtw.toHexString);
Tree.FullExpand;
Tree.Selected:=Tree.Items[0];
//Tree.SetFocus;
end;

{*******************************************************************************}
{constructor TBLS12CurvePairing.Create(inParametres:BLS12CurvesParamsDefinition);
begin
CurveParams:=nil;
SetStandardCurveParametres(inParametres);
CoordinatesSystem:=csProjective;
Fp12PoweringMode:=pmKarbina;
SetLoopMode(lpmAuto);
PairingAlgorithm:=bls12pOptAte;
end;      }

{*******************************************************************************}
constructor TBLS12CurvePairing.Create;
begin
inherited create(Aowner);
CurveParams:=nil;


SetStandardCurveParametres(sc128_1);
if CurveParams=nil  then exit;

_A:=Curveparams.A;
_B:=Curveparams.B;
_u:=Curveparams.u;
CoordinatesSystem:=csProjective;
Fp12PoweringMode:=pmKarbina;
SetLoopMode(lpmAuto);
PairingAlgorithm:=bls12pOptAte;
end;

{*******************************************************************************}
destructor TBLS12CurvePairing.Destroy;
begin
Dispose(CurveParams);
inherited;
end;

{*******************************************************************************}
function TBLS12CurvePairing.getA: String;
begin
Result:=_A.ToHexString;
end;

{*******************************************************************************}
function TBLS12CurvePairing.getAtw: String;
begin
Result:=CurveParams^.AtwFp4.toHexString;
end;

{*******************************************************************************}
function TBLS12CurvePairing.getB: String;
begin
Result:=_B.ToHexString;
end;

{*******************************************************************************}
function TBLS12CurvePairing.getBeta:integer;
begin
Result:=CurveParams^.FieldParam.Beta;
end;

{*******************************************************************************}
function TBLS12CurvePairing.getBtw: String;
begin
Result:=CurveParams^.Btw.toHexString;
end;

{*******************************************************************************}
function TBLS12CurvePairing.getCoord: CoordinatesSystem;
begin
Result:=CurveParams^.CoordSys;
end;

{*******************************************************************************}
function TBLS12CurvePairing.GetDefautG1Generator: FpPoint;
begin
Result.SetToDefaultGenerator;
end;

{*******************************************************************************}
function TBLS12CurvePairing.GetDefautG2Generator: Fp2Point;
begin
Result.SetToDefaultGenerator;
end;

{*******************************************************************************}
function TBLS12CurvePairing.getGamma: String;
begin
if  TwistMode=twDType then Result:=CurveParams^.TowerParam2.Gamma.toHexString
else Result:=CurveParams^.TowerParam2.Gamma.Inverse.toHexString
end;

{*******************************************************************************}
function TBLS12CurvePairing.getH: string;
begin
Result:=CurveParams^.H.toHexString;
end;

{*******************************************************************************}
function TBLS12CurvePairing.getHtw: string;
begin
Result:=CurveParams^.Htw.toHexString;
end;

{*******************************************************************************}
function TBLS12CurvePairing.getLp: string;
begin
Result:=CurveParams^.Lp.toHexString;
end;

{*******************************************************************************}
function TBLS12CurvePairing.getN: string;
begin
Result:=CurveParams^.N.toHexString;
end;

{*******************************************************************************}
function TBLS12CurvePairing.getP: string;
begin
Result:=CurveParams^.P.toHexString;
end;

{*******************************************************************************}
function TBLS12CurvePairing.getR: string;
begin
Result:=CurveParams^.R.toHexString;
end;

function TBLS12CurvePairing.GetRandomG1Point: FpPoint;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsRandomTorsionPoint;
end;

function TBLS12CurvePairing.GetRandomG2Point: Fp2Point;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsRandomTorsionPoint;
end;

{*******************************************************************************}
function TBLS12CurvePairing.getRtw: String;
begin
Result:=CurveParams^.Rtw.toHexString;
end;

{*******************************************************************************}
function TBLS12CurvePairing.getSigma: String;
begin
Result:=CurveParams^.TowerParam.Sigma.toHexString
end;

{*******************************************************************************}
function TBLS12CurvePairing.getTr: string;
begin
Result:=CurveParams^.Tr.toHexString;
end;

function TBLS12CurvePairing.GetTwm: TTwistModel;
begin
Result:=CurveParams.TwistMode;
end;

{*******************************************************************************}
function TBLS12CurvePairing.getu: String;
begin
Result:=_u.ToHexString;
end;

function TBLS12CurvePairing.HashToG1Point(id: string): FpPoint;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsTorsionFromString(id);
end;

function TBLS12CurvePairing.HashToG2Point(id: String): Fp2Point;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsTorsionFromString(id);
end;

{*******************************************************************************}
procedure TBLS12CurvePairing.RecomputeParametres;
var tmp:BLS12CurvesParamsDefinition;
begin
tmp.u:=_u.toHexString;
tmp.Beta:=beta;
tmp.Sigma:=Sigma;
tmp.Gamma:=Gamma;
tmp.A:=_A.toHexString;
tmp.B:=_B.toHexString;
tmp.SecurityLevel:=SecurityLevel;
ComputeBLS12Parametres(tmp,CurveParams);
end;

{*******************************************************************************}
procedure TBLS12CurvePairing.SetStandardCurveParametres(inParametres:StandardBLS12Curves);
begin
//ComputeBLS12Parametres(inParametres,CurveParams);
case inParametres of
sc128_1:FastLoadParams(sc128_1,CurveParams);
sc128_2_raz:FastLoadParams(sc128_2_raz,CurveParams);
sc128_3:FastLoadParams(sc128_3,CurveParams);
sc128_4:FastLoadParams(sc128_3,CurveParams);
end;

if CurveParams=nil then exit;

_A:=CurveParams.A;
_B:=curveparams.B;
_u:=curveparams.u;
if HammingWeight(Lp,true)<HammingWeight(lp,false) then PerferedPoweringMode:=lpmNaf
else PerferedPoweringMode:=lpmBinary;
end;

{*******************************************************************************}
procedure TBLS12CurvePairing.SetA(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TBLS12CurvePairing.SetB(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TBLS12CurvePairing.SetBeta(Value: Integer);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TBLS12CurvePairing.setCoord(value: CoordinatesSystem);
begin
CurveParams^.CoordSys:=value;
end;

{*******************************************************************************}
procedure TBLS12CurvePairing.SetCustomCurveParametres(inParametres: BLS12CurvesParamsDefinition);
begin
ComputeBLS12Parametres(inParametres,CurveParams);
end;

{*******************************************************************************}
procedure TBLS12CurvePairing.SetGamma(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TBLS12CurvePairing.SetSigma(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TBLS12CurvePairing.Setu(Value: String);
begin
RecomputeParametres;
end;

{********* Pairing function :Compute the pairing according to the specified algorithm***********}
function TBLS12CurvePairing.Paire(P: FpPoint; Q: Fp2Point): Fp12Int;
begin
if (P.Infinity)or(Q.Infinity) then raise Exception.Create('Can''t compute pairing for infinit points ....');
case PairingAlgorithm of
  bls12pOptAte:Result:=OptAtePairing(P,Q);
end;
end;


{******************** Set the loop mode of the pairing :Binary / Negatve Representation*************}
procedure TBLS12CurvePairing.SetLoopMode(Value: LoopPoweringMode);
begin
_LoopMode:=Value;
case value of
lpmBinary:begin
          case PairingAlgorithm of
          bls12pOptAte:Loop:=CurveParams.LoopBin;
          end;
          end;
lpmNaf:begin
       case PairingAlgorithm of
       bls12pOptAte:Loop:=CurveParams.LoopNaf;
       end;
      end;
lpmAuto:begin
     if CurveParams=nil then exit;
       case PairingAlgorithm of
        bls12pOptAte:begin
                 if ArrayHammingWeight(CurveParams.LoopNaf^)<ArrayHammingWeight(CurveParams.LoopBin^) then
                 Loop:=CurveParams.LoopNaf
                 else Loop:=CurveParams.LoopBin;
                 end;
       end;
       end;
end;
end;

procedure TBLS12CurvePairing.SetPArams(const Value: StandardBLS12Curves);
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

{******************** Final Exponentiation Step :Compute f^((p^12-1)/r) ***********************}
Procedure TBLS12CurvePairing.FinalPowerOptimalAteBLS12(ff:Fp12Int;u:LInt;PoweringMode:FpPoweringMode;var Result:Fp12Int); // Compute f^((p^12-1)/r)
var {xA,t0,t1,tmp,}f:Fp12Int;
    t:array[0..6]of Fp12Int;
    uu:lint;
begin
///***** Soft Exponentiation ***///
_Conjugate_FP12(ff,t[1]);
_Inv_FP12(ff,t[2]);
_Mul_FP12(t[1],t[2],f);
_Pow_FP12_P2(f,t[3]);
_Mul_FP12(t[3],f,f);
///***** Hard Exponentiation ***///
uu:=u-1;
_Sqr_FP12(f,t[0]);
_Mul_FP12(f,t[0],t[0]);
_Pow_FP12(f,uu,t[1],PoweringMode);
_Pow_FP12(t[1],uu,t[2],PoweringMode);
_Pow_FP12(t[2],u,t[3],PoweringMode);
_Mul_FP12(t[2],t[3],t[4]);
_Pow_FP12(t[4],uu,t[5],PoweringMode);
_Pow_FP12(t[5],u,t[6],PoweringMode);
_Mul_FP12(t[6],t[0],t[6]);
_Pow_FP12_P(t[5],t[5]);
_Mul_FP12(t[6],t[5],t[6]);
_Pow_FP12_P2(t[3],t[3]);
_Mul_FP12(t[6],t[3],t[6]);
_Pow_FP12_P3(t[2],t[2]);
_Mul_FP12(t[6],t[2],Result);
end;

{******************** Compute the Optimal Ate Pairing ***********************}
Function TBLS12CurvePairing.OptAtePairing(Pt:FpPoint;Qt:Fp2Point):Fp12Int;
var T,Q1,mQt,Q2,Q3,S:Fp2Point;
    i:integer;
    f:Fp12Int;
begin
f.SetTowerParams(CurveParams.TowerParam);
f.SetToOne;
T:=Qt;
_Neg_Fp2_Point(Qt,mQt);
T.SetPairingPointCoordinates(Pt.X,Pt.Y);
Q1.SetPairingPointCoordinates(Pt.X,Pt.Y);
Q2.SetPairingPointCoordinates(Pt.X,Pt.Y);
Q3.SetPairingPointCoordinates(Pt.X,Pt.Y);
S.SetPairingPointCoordinates(Pt.X,Pt.Y);
T.ComputeLigneValue:=true;
S.ComputeLigneValue:=true;
for i:=Length(Loop^)-2 Downto 0 do begin
                        case CoordinatesSystem of
                          csAffine:_Double_Affine_Fp2_Point(T,T,true);
                          csJacobian:_Double_Jacobian_Fp2_Point(T,T,True);
                          csProjective:_Double_Projective_Fp2_Point(T,T,True);
                        end;
                        _Sqr_FP12(f,f);
                        _Sparse_Mul_FP12(f,T.LineAtP,f,TwistMode);
                        if  Loop^[i]=1 then begin
                                         case CoordinatesSystem of
                                         csAffine:_Add_Affine_Fp2_Point(Qt,T,T,True);
                                         csJacobian:_Add_Jacobian_Fp2_Point(Qt,T,T,True);
                                         csProjective:_Add_Projective_Fp2_Point(Qt,T,T,True);
                                         end;
                                         _Sparse_Mul_FP12(f,T.LineAtP,f,TwistMode);
                                         end
                        else if (Loop^[i]=-1) then begin
                                                case CoordinatesSystem of
                                                csAffine:_Add_Affine_Fp2_Point(mQt,T,T,True);
                                                csJacobian:_Add_Jacobian_Fp2_Point(mQt,T,T,True);
                                                csProjective:_Add_Projective_Fp2_Point(mQt,T,T,True);
                                                end;
                                                _Sparse_Mul_FP12(f,T.LineAtP,f,TwistMode);
                                                end;
                        end;
FinalPowerOptimalAteBLS12(f,Qt.CurveParams.u,Fp12PoweringMode,Result);
end;

procedure Register;
begin
  RegisterComponents('Pairings', [TBLS12CurvePairing]);
end;

end.
