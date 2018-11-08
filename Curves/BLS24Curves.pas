unit BLS24Curves;

interface

uses  Vcl.ComCtrls, System.SysUtils,VLargeIntegers,LargeIntegers,Fp12Arithmetic,Fp6Arithmetic,Fp4Arithmetic,Fp2Arithmetic
      ,GeneralTypes,VCL.dialogs,System.classes,ECCFp4,ECCFp,Fp24Arithmetic,TreeView1;


Type
   GtBls24=Fp24Int;
   G1BLS24=FpPoint;
   G2BLS24=Fp4Point;
   LargeInt=Lint;
   BLS24CurvesPairingAlgos=(blspOptAte);
   StandardBLS24Curves=(scBLS24at256_0,scBLS24at256_1,scBLS24at256_2,scBLS24Razat256,scBLS24at192_1,scBLS24at192_2,scBLS24at192_3,scBLS24at192_Raz,scBLS24at320);

   BLS24CurvesParamsDefinition=record
                            SecurityLevel:String;
                            u:String;  // the paramater of generation for the BLS24 curve
                            Beta:Integer; // non-square elements of the irredictible polynomial on FP2
                            Sigma:string; // non-square non-cube elements of the irredictible polynomial on FP4
                            Gamma:string; // non-square non-cube elements of the irredictible polynomial on FP8
                            A,B:String; // parametres of the curve   y^2=x^3+Ax+B
                            TwistMode:TTwistModel;
                            GenratorX:String;
                            TwistGeneratorSeed:Word;
                            end;
   ListOfBLS24Params=array of BLS24CurvesParamsDefinition;

 TBLS24CurvePairing=class (TCurve)
                      {Definition of a Curve with respect the the Weistrass model   Y^2=X^3+A*X+B (mod P)}
              private
              _A,_B:Lint;
              _u:Lint;
              _LoopMode:LoopPoweringMode;                     //  Loopmode :Negative/Binary Representation
              _Fp24PoweringMode:FpPoweringMode;             // Final Exponentiation Powering Mode : Beuchat Approach, Karabina Approach, Classical Square and Multiply
              _PairingAlgorithm:BLS24CurvesPairingAlgos;         //  Pairing Algorithme OptAte,R-Ate.....
              params:StandardBLS24Curves;
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
              Procedure FinalPowerOptimalAteBLS1(f:Fp24Int;u:LInt;PoweringMode:FpPoweringMode;var Result:Fp24Int); // Compute f^((p^24-1)/r)
              Function OptAtePairing(Pt:FpPoint;Qt:Fp4Point):Fp24Int;  // Compute Optimal Ate Pairings
              procedure SetLoopMode(Value:LoopPoweringMode);
              procedure SetPArams(const Value: StandardBLS24Curves);
            public
              SecurityLevel:String;    /// A symbolic identifier of the curve
              PerferedPoweringMode:LoopPoweringMode;
              Loop:PLIntArrayForm;
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
              function Paire(P:FpPoint;Q:Fp4Point):Fp24Int;
              Constructor Create(AOwner : TComponent); override;
              destructor Destroy;
              procedure SetStandardCurveParametres(inParametres:StandardBLS24Curves);
              procedure SetCustomCurveParametres(inParametres:BLS24CurvesParamsDefinition);
              procedure GenerateParamsTreeView(Tree:TTreeView);override;
              function GetRandomG1Point:FpPoint;
              function GetRandomG2Point:Fp4Point;
              function GetDefautG1Generator:FpPoint;
              function GetDefautG2Generator:Fp4Point;
              function HashToG1Point(id:string):FpPoint;
              function HashToG2Point(id:String):Fp4Point;
              published
              property Fp24PoweringMode:FpPoweringMode read _Fp24PoweringMode write _Fp24PoweringMode;
              property PairingAlgorithm:BLS24CurvesPairingAlgos read _PairingAlgorithm write _PairingAlgorithm ;
              property LoopMode:LoopPoweringMode read _LoopMode write SetLoopMode;
              property CoordinatesSystem:CoordinatesSystem read getCoord write SetCoord;
              property Parametres:StandardBLS24Curves read params write Setparams;
            end;


const

  BLS24_256_0_Params:BLS24CurvesParamsDefinition=(SecurityLevel:'256bit';u:'0x7FFF804000000000';Beta:-1;sigma:'1+u*1';Gamma:'(0+u*0,1+u*0)';A:'0';B:'4';TwistMode:twMType;GenratorX:'5';TwistGeneratorSeed:55);
  BLS24_256_1_Params:BLS24CurvesParamsDefinition=(SecurityLevel:'256bit';u:'0xE000000000058400';Beta:-1;sigma:'1+u*1';Gamma:'(0+u*0,1+u*0)';A:'0';B:'-2';TwistMode:twDType;GenratorX:'3';TwistGeneratorSeed:55);
  BLS24_256_2_Params:BLS24CurvesParamsDefinition=(SecurityLevel:'256bit';u:'-0x10000140000800000';Beta:-1;sigma:'1+u*1';Gamma:'(0+u*0,1+u*0)';A:'0';B:'-2';TwistMode:twDType;GenratorX:'3';TwistGeneratorSeed:55);
  BLS24_Razvan_256_Params:BLS24CurvesParamsDefinition=(SecurityLevel:'256bit';u:'-0x9FFFFFFFEFFFFC000000000000';Beta:-1;sigma:'1+u*1';Gamma:'(0+u*0,1+u*0)';A:'0';B:'-2';TwistMode:twDType;GenratorX:'3';TwistGeneratorSeed:55);
  //  -1+2^19-2^24+2^27-2^48 (Proposed By Craig Castello) Hamming Weight=5
  BLS24_192_1_Params:BLS24CurvesParamsDefinition=(SecurityLevel:'192bit';u:'-0xFFFFF8F80001';Beta:-1;sigma:'1+u*1';Gamma:'(0+u*0,1+u*0)';A:'0';B:'1';TwistMode:twDType;GenratorX:'6';TwistGeneratorSeed:55);
  //  2^31-2^33-2^49 (Proposed By Craig Castello) Hamming Weight=3
  BLS24_192_2_Params:BLS24CurvesParamsDefinition=(SecurityLevel:'192bit';u:'0x3FFFF7FF0000';Beta:-1;sigma:'1+u*1';Gamma:'(0+u*0,1+u*0)';A:'0';B:'-2';TwistMode:twDType;GenratorX:'3';TwistGeneratorSeed:55);
  //   Proposed by Loubna Ghammam , Hamming weight=3 but Sigma is not suitable for fast computation
  BLS24_192_3_Params:BLS24CurvesParamsDefinition=(SecurityLevel:'192bit';u:'0xFFFFC4000000';Beta:-1;sigma:'3+u*1';Gamma:'(0+u*0,1+u*0)';A:'0';B:'5';TwistMode:twDType;GenratorX:'-1';TwistGeneratorSeed:55);
  BLS24_Razvan_192_Params:BLS24CurvesParamsDefinition=(SecurityLevel:'192bit';u:'-0x10007FFFFFFFE40';Beta:-1;sigma:'1+u*1';Gamma:'(0+u*0,1+u*0)';A:'0';B:'-2';TwistMode:twDType;GenratorX:'-1';TwistGeneratorSeed:55);
  BLS24_320_1_Params:BLS24CurvesParamsDefinition=(SecurityLevel:'320bit';u:'0x7FFFFFFFFFFFFFFFFC000024000';Beta:-1;sigma:'1+u*1';Gamma:'(0+u*0,1+u*0)';A:'0';B:'-2';TwistMode:twDType;GenratorX:'3';TwistGeneratorSeed:55);



  procedure ComputeBLS24Parametres(Params: BLS24CurvesParamsDefinition; var Result:PtrCurveParams);
  procedure Register;

implementation



Procedure FastLoadParams(param:StandardBLS24Curves;var Result:PtrCurveParams);
var rho,tmp,e:Lint;
    i:integer;
    TwoinFp2:Fp2Int;
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
scBLS24at256_0:begin
               SecurityLevel:='256bit';
               Family:=cfBLS24;
               u:='0x7FFF804000000000';
               P:='0x1554806E66E1A9B4339DBF675F5F3308BC8E6E6BB24AB1FC1E7F3CDFBB8356CA3581F90268DAFD557A24AE4684807DABE7229A7F26B6358AB2A2AD5500000015552AC0154005557FFFD56AAAAAAAAB';
               new(p.InverseFordivision);
               p.Limit:=2*p.Data.i16[-2];
               p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
               p.id:='Y';
               r:='0xFFF8041BE3CF541B7B316894524F57B159B5EE437E506FF000FFFFFFFFFFFFF0003FDFA06027A02017E807FF00000000000000000000000000000000000001';
               Tr:='0x7FFF804000000001';
               Lp:= '0x7FFF804000000000';
               n:='0x1554806E66E1A9B4339DBF675F5F3308BC8E6E6BB24AB1FC1E7F3CDFBB8356CA3581F90268DAFD557A24AE4684807DABE7229A7F26B6358AB2A2AD5500000015552AC0154005550000552AAAAAAAAB';
               htw:='0x328B135818BC5B0498C75D1FEA2C5DC2F21383A8EB9F1BC494738D191896D031117FDE40BA621396B39C756C70E5B5C3C6B932EC51070B35C009CAB7B5999821F888635562C110DAA52B196E5733C114986F8EDD45D991F295EB6950845C'+
               'FB875E4BE1D93E6E85559805405BA6E29CAA1CFA54B2A01510B1BE96EED703A4C03D6E4BCC3255829AA059573A4CCE98F88D606B0BD14BAB305ED44E430EF2A08FCAF9EE5FD2A5BF273FD1300A50C37A9326144824D844ACA55BCB2A175273DB7745B'+'DFB02142D41C0DFFC96108B431DDA9C098CDFF17A2F641FD0223EEB76BDA3CE0024084F0E8703CC83C120AFD9EB0D9F0554CA497B3E6D6E9E06524';
               rtw:=r;
               A:=0;
               B:=4;
               H:='0x7FFF803FFFFFFFFF';
               FieldParam:=InitFieldParams(-1,p);
               New(FieldParam.MontgomeryData);
               InitMontgomeryStruct(P,FieldParam.MontgomeryData);
               New(TowerParam);
               TowerParam^.FieldParam:=FieldParam;
               TowerParam^.Sigma.SetFieldParams(TowerParam^.FieldParam);
               TowerParam.Sigma.SetFormString('1+u*1');
               TowerParam.pmod8:=p mod 8;
               TwoinFp2.SetFieldParams(TowerParam^.FieldParam);
               TwoinFp2.SetFormString('2+u*0');       // for squaring over FP4
               TowerParam^.inv_2.SetFieldParams(TowerParam^.FieldParam);
               TowerParam^.inv_2:=TwoinFP2.Inverse;
               New(TowerParam2);
               TowerParam2.Tower8paprms:=TowerParam;
               TowerParam2.pmod8:=TowerParam.pmod8;
               TowerParam2.Sigma:=Towerparam.Sigma;
               TowerParam2.FieldParam:=FieldParam;
               TowerParam2.Gamma.SetTowerParams(TowerParam);
               TowerParam2.Gamma.SetFromStrings('0+u*0','1+u*0');
               BtwFp4.SetTowerParams (TowerParam);
               TwistMode:=twMType;
               BtwFp4.a.SetFormString('0+u*0');
               BtwFp4.b.SetFormString('4+u*0');
               LoopBin^:=Lp.ToIntArray;
               LoopNaf^:=Lp.ToNafArray;
               LoopTateBin^:=n.ToIntArray;
               LoopTateNaf^:=n.ToNafArray;
               BasePointY:='0x7C5FC9C51B41750B7B3F17CF2BD07EE2EF8692DA3529D67A6CBE59ACD5D5E95B0827293816F76EE3579504EA971067F55294AB9CA52581A80009D5433F34CEFF55DA156D5073AC79F6DC6FCFA1B44';
               BasePointX:='0x3';
               BLS24TwistGeneratorX.a.SetFormString('0xECAAE17829214D7E4AC348930470CCFF5F12F6CC9294A546248503EBAA0820EBDC88EAD05DCEA28F5E34A79086EE052BF55B4D32625F7432EFE2C38D03E7FE3E0FF0BFA31CDEFB534EC88B1ECB84D'+'+u*0x59DC35FEBD95B0FBA65C387330B4C636BDA2098295EAF2D35D2BE65E6961F02117193FE13ED6D67A9FAF840BCFD893EF9795E3D171281536E9EA352850D4195852C009A769A158815B342A279866F');
               BLS24TwistGeneratorX.b.SetFormString('0x13F1ADA512DCF9A576761110788B90144EF6A65A72DD5A7E5A20FB6491E3BA22DE90331CA3253F2FA4DF71B12A8963654A3FDB683DD0DE2E530DF968ACB5292D3EA4CC086A8E5DC1C588F4FC9773A'+'+u*0x15205B18E0E089F84A28E817AEEA0AF253699D7BA55DCD43D1CF54CA173133B1BC8E1B485F114AC2C73288228753B372165455FFA60877C49895A1821B564BF742A70C9840F99ED7CB14FD2A0DB89B');
               BLS24TwistGeneratorY.a.SetFormString('0x18D15AFD8686CD838F993964A357406C7A0CDAE56A315AC2BAB3EE3D94C934FB0FCF2875942836CA1C57CBBCA854E88AD25E6594ADF7D65B6CD8F6C1D64B3556A8E87940C46949DEACB72DFCCA6E0'+'+u*0x8427D3D4CEE526197EBE66849D1B3762752FEE5665FA0F2C1B3A9130C7BCFC769019F2FE0A4A2989E0006CA677C3A01A6116F92D6BB59D716497CFFED1BA9B15B465AFB92DE5682539F979C4B910D');
               BLS24TwistGeneratorY.b.SetFormString('0x9715DFF2649109C6B54D93F577AFF13F9E01207D7D7EC62F13ADC4BAA77E1CEC885694E4054C217AC877D0FF211B7A040DB76145FF06C448744B8124CBBA6696CC44AE6E97C435A4F9B280F7DF00C'+'+u*0x122274B1EA2F5DB32C6F0D044D6DF3CBC6317407CC7649650B04459D9C731B8FE5B7E36AA2A50DF4E35A9B48B892778EAA70F65459F36B9DFF9102A96FE959CF4FD23BEBA97F1A62C3E43CB16993EE');
        end;
scBLS24at256_1:begin
               SecurityLevel:='256bit';
               Family:=cfBLS24;
               u:='0xE000000000058400';
               P:='0x1672F941555ADC6DB8DBBAB222012CF069DB4EDCD279A04B8128522FC8536700E905C47818A293417DBD0ED4C9E9BDEB04944B45ED24C4AD6169FEA27B9CA69426C706A9CC83D2BA41E37FB4CEB1D6AB';
               new(p.InverseFordivision);
               p.Limit:=2*p.Data.i16[-2];
               p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
               p.id:='Y';
               r:='0x57F6C100001153E49700017E4AE3437012D39DD5B6C49457FF15181B4C13865780BD45C6AE3924ABAD6721D3D30F6BB01CCA4E0C1B00FC627EA27F0000000001';
               Tr:='0xE000000000058401';
               Lp:= '0xE000000000058400';
               n:='0x1672F941555ADC6DB8DBBAB222012CF069DB4EDCD279A04B8128522FC8536700E905C47818A293417DBD0ED4C9E9BDEB04944B45ED24C4AD6169FEA27B9CA69426C706A9CC83D2B961E37FB4CEAC52AB';
               htw:='0xB474DEDA1F68ECC4C6A003DFB5B846962C26FB19FDC8F36C8C7FD184773A15F791C6195C5C6FE5DC6A522223CFBCFC6D80091B6D6F8651052C095660DD29091156997619DF0B942C9C1CC2996B0500E'
               +'F943FC0E376B9E97FDA595FF9B4F8FDC321331D7BDD7511C14CA7B8FA8E30189683903A3462389EA8EBC969305DDFB571FFFA898246954F1C77F6CCEB73762B35DE775150C7B1A127659560C4650FAEECC4E2'
               +'A3EF131622380F9240A7CFCBAB42C957F9AE28F0E3C8EC5D488A6E45D1960694020985DBC4D7EDE71B7AA2E6F220F8A033503C4778FCA71844B79062001D569D3ED4426CAC1BB2DBD97171E9A997A5EA6D561'
               +'772C15A3E8ADB8AB0524';
               rtw:=r;
               A:=0;
               B:=-2;
               H:='0xE0000000000583FF';
               FieldParam:=InitFieldParams(-1,p);
               New(FieldParam.MontgomeryData);
               InitMontgomeryStruct(P,FieldParam.MontgomeryData);
               New(TowerParam);
               TowerParam^.FieldParam:=FieldParam;
               TowerParam^.Sigma.SetFieldParams(TowerParam^.FieldParam);
               TowerParam.Sigma.SetFormString('1+u*1');
               TowerParam.pmod8:=p mod 8;
               TwoinFp2.SetFieldParams(TowerParam^.FieldParam);
               TwoinFp2.SetFormString('2+u*0');       // for squaring over FP4
               TowerParam^.inv_2.SetFieldParams(TowerParam^.FieldParam);
               TowerParam^.inv_2:=TwoinFP2.Inverse;
               New(TowerParam2);
               TowerParam2.Tower8paprms:=TowerParam;
               TowerParam2.pmod8:=TowerParam.pmod8;
               TowerParam2.Sigma:=Towerparam.Sigma;
               TowerParam2.FieldParam:=FieldParam;
               TowerParam2.Gamma.SetTowerParams(TowerParam);
               TowerParam2.Gamma.SetFromStrings('0+u*0','1+u*0');
               BtwFp4.SetTowerParams (TowerParam);
               TwistMode:=twDType;
               BtwFp4.a.SetFormString('0+u*0');
               BtwFp4.b.SetFormString('0x1672F941555ADC6DB8DBBAB222012CF069DB4EDCD279A04B8128522FC8536700E905C47818A293417DBD0ED4C9E9BDEB04944B45ED24C4AD6169FEA27B9CA69426C706A9CC83D2BA41E37FB4CEB1D6AA+u*1');
               LoopBin^:=Lp.ToIntArray;
               LoopNaf^:=Lp.ToNafArray;
               LoopTateBin^:=n.ToIntArray;
               LoopTateNaf^:=n.ToNafArray;
               BasePointY:='5';
               BasePointX:='3';
               BLS24TwistGeneratorX.a.SetFormString('0x699DC334F81547AFEB4FA9E1DFFCD39014E2066ABC972BDB7688151F504781D38EA4F60DFA6F4E63AF7711AA496537D8B3FD0F85AD967744D25544ED67637849033DE9ADDAAEA9AA03FBC8DC511817A'+'+u*0x105E6AC86B3063FE0202E630CE53213438B60249B25C5724CDB98B188F58541BC366A64F984352F1F89F68DF4B57DE499457919BDB6E202A0D576BA08AC867F0D7E67CD33B54FF2F1175A114080FD98A');
               BLS24TwistGeneratorX.b.SetFormString('0xD37423ABF2CF11F2586AB2C126712AE58C238583945AAD303F5398669240032AFEA3EE59F04B772B015D99B3B71CBEF31AD686F3236D6F62C3374DD4AFC54C9B2C55A2224D8726007C40C75356444F9'+'+u*0x47563630EBFB1570D869355E84EFE0F902AE9B04F8AB03C21490438A677D1D817D37DC9B76BE66DD1186000BEEC47D4297BC5CBBA70DF7E6BC39BE7A918A8BB1A831A11E3A895DB7D75AC530418ECE1');
               BLS24TwistGeneratorY.a.SetFormString('0x3ED76EE930A340938A699DE06DD4B5A72291C3C27897271D285F0F68D6096F1E2AC9A6C608CC9E36176DD2D57E754A41C13DE879BF25FAF6D62C9CC97141C621D659BEB7A405481AB03382283C975BE'+'+u*0x660150D3A5600DC8C8BDE3A5B3E57D15E7E1C881FE1ED547A54036069A5AEACF267CB50700566880A2C179B639003EB0B0556425A4E627D997A8915E62ED55E860EBAB57799930843DB1F67EC3009A9');
               BLS24TwistGeneratorY.b.SetFormString('0xEEFF0416954D10E11614907D39746257B567204EA8EF5AA82C18B2901865B16F7F996D5F8BA31AFB7422EB4372CB077A857F39B7FBDC6F5A99DFE402BE5EEC2F00F65B6DD9783877D8EF59F3011A98E'+'+u*0x38B9C100E6C5A0C3B2B96C9A26131FED606D49475CC82D47030E57B848705C5B92F7700CD852B5C95AE0159019D7252773A7A81AA762617F1EB41089BFE2707F4EBAAD1C775EEE2F2E4869ED5614A1D');
        end;
scBLS24at256_2:begin
               SecurityLevel:='256bit';
               Family:=cfBLS24;
               u:='-0x10000140000800000';
               P:='0x55559800191AB0B956623378F26912D55C5963C0A540368AB35037113960F82A8188A09F59D148A8E21765C27A04F2C2AEA20EC46213371754C3CFA7BF2AA438554CB7FFFB85555BAAAAB955552AAAAB';
               new(p.InverseFordivision);
               p.Limit:=2*p.Data.i16[-2];
               p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
               p.id:='Y';
               r:='0x10000A0002FC009060135261F0D626BABE63AEBE2DEB41D759AEAF3E1AC13525F482FB0BEF4A13F0B00F2AEFF857FFD27FFF5FFFFF00000000000000000000001';
               Tr:='-0x100001400007FFFFF';
               Lp:= '0x10000140000800000';
               n:='0x55559800191AB0B956623378F26912D55C5963C0A540368AB35037113960F82A8188A09F59D148A8E21765C27A04F2C2AEA20EC46213371754C3CFA7BF2AA438554CB7FFFB85555CAAAACD5555AAAAAB';
               htw:='0x3291E065BF03A858CB3E09CA09D9C7D06B3D379BD319AC4DB25C18FDE2555D5CD465B028F6E3716E4C1C5DCD28504D9B0F23B9C1A1E07FE1FBFCA8E70BC433D79E937B0FC3A2108055A70F7A8A8BB7AD'
               +'4E30CE419CDE6F4BBCA549452CCBF28EFA84F3FBD58E2A134E6A855CEE3BED72AADD871AE8F0AF71230ABC00172FE2697E59B1D8541B6382697529B232BEC807DAB8FDA2BA3677EB797AFA05CC4316F91BEDB3'
               +'0295EC253BA7B03DF0DCFB013388E49D19BEEE7F541EDCD4C0DF3187AFB240D175CD1EC05ADD19407C745C1D5BF8ADC8A2F498CACE62DE88BD1DD28EEAE598A1A9691022CE831859FA5736BB084D4F7B43CBB4269948AE32915E06524';
               rtw:=r;
               A:=0;
               B:=-2;
               H:='-0x10000140000800001';
               FieldParam:=InitFieldParams(-1,p);
               New(FieldParam.MontgomeryData);
               InitMontgomeryStruct(P,FieldParam.MontgomeryData);
               New(TowerParam);
               TowerParam^.FieldParam:=FieldParam;
               TowerParam^.Sigma.SetFieldParams(TowerParam^.FieldParam);
               TowerParam.Sigma.SetFormString('1+u*1');
               TowerParam.pmod8:=p mod 8;
               TwoinFp2.SetFieldParams(TowerParam^.FieldParam);
               TwoinFp2.SetFormString('2+u*0');       // for squaring over FP4
               TowerParam^.inv_2.SetFieldParams(TowerParam^.FieldParam);
               TowerParam^.inv_2:=TwoinFP2.Inverse;
               New(TowerParam2);
               TowerParam2.Tower8paprms:=TowerParam;
               TowerParam2.pmod8:=TowerParam.pmod8;
               TowerParam2.Sigma:=Towerparam.Sigma;
               TowerParam2.FieldParam:=FieldParam;
               TowerParam2.Gamma.SetTowerParams(TowerParam);
               TowerParam2.Gamma.SetFromStrings('0+u*0','1+u*0');
               BtwFp4.SetTowerParams (TowerParam);
               TwistMode:=twDType;
               BtwFp4.a.SetFormString('0+u*0');
               BtwFp4.b.SetFormString('0x55559800191AB0B956623378F26912D55C5963C0A540368AB35037113960F82A8188A09F59D148A8E21765C27A04F2C2AEA20EC46213371754C3CFA7BF2AA438554CB7FFFB85555BAAAAB955552AAAAA+u*1');
               LoopBin^:=Lp.ToIntArray;
               LoopNaf^:=Lp.ToNafArray;
               LoopTateBin^:=n.ToIntArray;
               LoopTateNaf^:=n.ToNafArray;
               BasePointY:='0x7C5FC9C51B41750B7B3F17CF2BD07EE2EF8692DA3529D67A6CBE59ACD5D5E95B0827293816F76EE3579504EA971067F55294AB9CA52581A80009D5433F34CEFF55DA156D5073AC79F6DC6FCFA1B44';
               BasePointX:='3';
               BLS24TwistGeneratorX.a.SetFormString('0x2C0A550DB093798C803A710523AAD2C4FB4755AA94E6848B56C45A385AEF9FDEFBC1AFA88BD9BB27367B827FFD6D6CBE14654EF2231C8FE4FD2F285EB5D950D6FC81D2654937D57B26006802C7BC0421'+'+u*0x36EA4A6B581891FAE5AE23AB2B53A7A38D3CE30F05AAB2714EF37CC2A936FC6CA708A9EE14F8B18AC2F0B02E03A139FD3CF5E430017AE1190B0015C4D4F53A3362E47A8841B8D01F77BD0E7F7AC5D772');
               BLS24TwistGeneratorX.b.SetFormString('0x2385294CFB02BBDAD77EE9FFF1AE9E0231C3FEAF740567B99940D94387F29AD11CD0B5418643E3813DAA7B4EE7EFB6015DBC3FA51F7768BC4010FEDB382B00B464BAF20EFC0C5E4B61F5ADD4519D5B1B'+'+u*0x497E34B9CB9DC48DA721BC4C7F27948060F25EFF17CB1195D7B011EC25A29565841A891D626E05F7722201B910A20E317874E7FB9A729F77C340C4F4BA27DEF1F8BB40B0081EC3991553563B19E48D18');
               BLS24TwistGeneratorY.a.SetFormString('0x4E19B74382AF3439C7E0F9675D66F33377E7F3D40F11F6CA28CACD8349135C1F3F2DB9EE2E7EE8846C160B8E350F6AD353F20CE9AB9D4AF4B76A4EEAD6C71449D96964B41A95D508E1F2B0F471A721E5'+'+u*0x531496E4AC5B69233E78D1D3EEA4ABA5E1ECF4068315CA871368202B0A98B168573F91C5705631576F70382DC571307D6787ADB66CE142BCC014F71CB8BA50FB3820852E6CF78D7C25D7A78E647E68FF');
               BLS24TwistGeneratorY.b.SetFormString('0x50A56A89AE66F5E8A521BA6A0A0C90CEE5B890CD072C2FE9C9D01B7E58AF6DA83E2E0A20EEC2F09AAB7AA5BE0592AD8CEAA2C821379D38D0DA05191661900EAEEA5685FD9B89496E85F873208EB57C61'+'+u*0x12672BAF52B1AE1E1D7F90016483173B9285A010235AFE5A8D3F59F65EA39EC9795F08EFA4C310E103A07F105B9D98E53652E3D56F4A6D61E7C7B1BE16754C0E6199C907A5FB06A30BF34B72130E53C6');
        end;
scBLS24at192_1:begin
               SecurityLevel:='192bit';
               Family:=cfBLS24;
               u:='-0xFFFFF8F80001';
               P:='0x55553DE5583EE8DE07A5BD8D1F2B1D181350B13B2AAAC13EA6DB72C4CB7C8EAF4D23C01812428A6B82580149B426A18F5C662A623EC8036A82D2AAAB';
               new(p.InverseFordivision);
               p.Limit:=2*p.Data.i16[-2];
               p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
               p.id:='Y';
               r:='0xFFFFC7C0057046B26BDC4CE1022439D3586C9E35253FDC6BB4C486B2C06E624FA8DCAE61A5547EAE9C43A57FE3E00001';
               Tr:='-0xFFFFF8F80000';
               Lp:= '0xFFFFF8F80001';
               n:='0x55553DE5583EE8DE07A5BD8D1F2B1D181350B13B2AAAC13EA6DB72C4CB7C8EAF4D23C01812428A6B82580149B426A18F5C662A623EC9036A7BCAAAAC';
               htw:='0x32913587F9AAF0815EC01514E97F5171A4CA912B1E92796165045C1EE46E29725C36BA5730A9966D134672A0ACD83DE241FC225B9A4DE4238BC30AF1A442F7D1A302'
               +'E778517DD4E0A5D8C0D28B557947D57CCA1193F79A6E3B996BA4CBA3D5BDC8630CA422A449625EC94B42B91C212AD05A020260B0BDD8CA1942611289DEA2103E0DCC354D1C'
               +'C6E668375BB5FAE7027B000265C179CD64BA230C4573D849156D1A36F224CE1A9735BC878AA76791448491CCB9CB8F2B5DC258BE55575BA79';
               rtw:=r;
               A:=0;
               B:=1;
               H:='-0xFFFFF8F80002';
               FieldParam:=InitFieldParams(-1,p);
               New(FieldParam.MontgomeryData);
               InitMontgomeryStruct(P,FieldParam.MontgomeryData);
               New(TowerParam);
               TowerParam^.FieldParam:=FieldParam;
               TowerParam^.Sigma.SetFieldParams(TowerParam^.FieldParam);
               TowerParam.Sigma.SetFormString('1+u*1');
               TowerParam.pmod8:=p mod 8;
               TwoinFp2.SetFieldParams(TowerParam^.FieldParam);
               TwoinFp2.SetFormString('2+u*0');       // for squaring over FP4
               TowerParam^.inv_2.SetFieldParams(TowerParam^.FieldParam);
               TowerParam^.inv_2:=TwoinFP2.Inverse;
               New(TowerParam2);
               TowerParam2.Tower8paprms:=TowerParam;
               TowerParam2.pmod8:=TowerParam.pmod8;
               TowerParam2.Sigma:=Towerparam.Sigma;
               TowerParam2.FieldParam:=FieldParam;
               TowerParam2.Gamma.SetTowerParams(TowerParam);
               TowerParam2.Gamma.SetFromStrings('0+u*0','1+u*0');
               BtwFp4.SetTowerParams (TowerParam);
               TwistMode:=twDType;
               BtwFp4.a.SetFormString('0+u*0');
               BtwFp4.b.SetFormString('0x2AAA9EF2AC1F746F03D2DEC68F958E8C09A8589D9555609F536DB96265BE4757A691E00C09214535C12C00A4DA1350C7AE3315311F6401B541695556+u*0x2AAA9EF2AC1F746F03D2DEC68F958E8C09A8589D9555609F536DB96265BE4757A691E00C09214535C12C00A4DA1350C7AE3315311F6401B541695555');
               LoopBin^:=Lp.ToIntArray;
               LoopNaf^:=Lp.ToNafArray;
               LoopTateBin^:=n.ToIntArray;
               LoopTateNaf^:=n.ToNafArray;
               BasePointY:='0x7C5FC9C51B41750B7B3F17CF2BD07EE2EF8692DA3529D67A6CBE59ACD5D5E95B0827293816F76EE3579504EA971067F55294AB9CA52581A80009D5433F34CEFF55DA156D5073AC79F6DC6FCFA1B44';
               BasePointX:='6';
               BLS24TwistGeneratorX.a.SetFormString('0x458C152D5BDC788F000E92D37B31F7BE74CE2FC511434CC22786F46A7A8552415300B74E1D0E5C049300E233E8972A313019C139C5C9F505A32400CE+u*0x38A299BBB92F488CACF54EAA0E654212682677B40F857E76FADA6BEAC0B0F619820E572C0050CC84F27512A466D61521362E6CDF3FC862A95CEC96');
               BLS24TwistGeneratorX.b.SetFormString('0x106E2757B192A9FE92E17638055CD78BD09F5ABF2CFA75DCBC7E7E4A6B97D80D9035FDF432177FE5357C4175A705ED9C4D5CD300C7C969F690929905+u*0x2B6193DBFF8C771B4C22DB25F8F239C5CDC1F848867760BB681BFA33DA6BC2D9773A95FE51EBAEDAF7796FDDF2B6D9751C4E572DF54696055F463149');
               BLS24TwistGeneratorY.a.SetFormString('0x310CE549AC63CC4220619E6F0E41A9214260FBF163F4B67F428EB10500FB068F0699A56450FF79121AEB882B2288DB9FBC1ACAEDB27FB04E38D883F0+u*0x4DEF1CEA9FC9A35D25AFF847D9C0417252CB6E767CC7E36EAB07CA3E2A6C77A2395B73FB0BF1E3519B36B2CB956B2FC469CDEC0526B0E4157977D226');
               BLS24TwistGeneratorY.b.SetFormString('0x2242A358AB7BA4B27C9F9B95FDE75D644887FFC2BF9BB27F8A0079414961A5A05395D00A4059481BC02AE522CA9CA4295E2FA33DC7E652002299A813+u*0x504F3B18A85613CF66F113CAC17504D81CCA6BE70539783466CB5CB61DC9F7325BCA55580F99090F97D45C0AFEF335DF6C4A2B1CE12BC479A9A51DF4');
        end;
scBLS24at192_2:begin
               SecurityLevel:='192bit';
               Family:=cfBLS24;
               u:='0x3FFFF7FF0000';
               P:='0x5554EA9D916197962E7B501109DC2A12B3A8FF044390BED40FA77ADDB490D3B64A4F077540840941833EEF370FD3F7F4954D56AC555A7FFAAAB';
               new(p.InverseFordivision);
               p.Limit:=2*p.Data.i16[-2];
               p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
               p.id:='Y';
               r:='0xFFFEFFE0701BE5B5831221F8C2B5C911410D128CF3E5D083F0EE16E83A02A0C007F8FE7FDFFF0000000000000001';
               Tr:='0x3FFFF7FF0001';
               Lp:='0x3FFFF7FF0000';
               n:='0x5554EA9D916197962E7B501109DC2A12B3A8FF044390BED40FA77ADDB490D3B64A4F077540840941833EEF370FD3F7F4954D56A8555B000AAAB';
               htw:='0x3290979C657C871698ED349F4DA284A2080995148B19F860AF5BE8DCF0EF418B1F22650E55CC7363EA9DE729B9B1DA3E0C3A6CB791F491F99C'
               +'ABAD041E18EB0E6B9B8DB9E18D06F76DA53966BBFC244288ABACA46BFBC3ADD89A4A5AD77399CFB48441B2E9276B9EA139E553305097E066E534EA7D'
               +'0B61FD740558EDE0CCBED3F4536DF23C0AE14278EF804D8DF7E96FC68A358EEB907D5E740C25BC18BAFCB4A7CD097C1935247B2C1053B3F54BC6003D15F9A4E386524';
               rtw:=r;
               A:=0;
               B:=-2;
               H:='0x3FFFF7FEFFFF';
               FieldParam:=InitFieldParams(-1,p);
               New(FieldParam.MontgomeryData);
               InitMontgomeryStruct(P,FieldParam.MontgomeryData);
               New(TowerParam);
               TowerParam^.FieldParam:=FieldParam;
               TowerParam^.Sigma.SetFieldParams(TowerParam^.FieldParam);
               TowerParam.Sigma.SetFormString('1+u*1');
               TowerParam.pmod8:=p mod 8;
               TwoinFp2.SetFieldParams(TowerParam^.FieldParam);
               TwoinFp2.SetFormString('2+u*0');       // for squaring over FP4
               TowerParam^.inv_2.SetFieldParams(TowerParam^.FieldParam);
               TowerParam^.inv_2:=TwoinFP2.Inverse;
               New(TowerParam2);
               TowerParam2.Tower8paprms:=TowerParam;
               TowerParam2.pmod8:=TowerParam.pmod8;
               TowerParam2.Sigma:=Towerparam.Sigma;
               TowerParam2.FieldParam:=FieldParam;
               TowerParam2.Gamma.SetTowerParams(TowerParam);
               TowerParam2.Gamma.SetFromStrings('0+u*0','1+u*0');
               BtwFp4.SetTowerParams (TowerParam);
               TwistMode:=twDType;
               BtwFp4.a.SetFormString('0+u*0');
               BtwFp4.b.SetFormString('0x5554EA9D916197962E7B501109DC2A12B3A8FF044390BED40FA77ADDB490D3B64A4F077540840941833EEF370FD3F7F4954D56AC555A7FFAAAA+u*1');
               LoopBin^:=Lp.ToIntArray;
               LoopNaf^:=Lp.ToNafArray;
               LoopTateBin^:=n.ToIntArray;
               LoopTateNaf^:=n.ToNafArray;
               BasePointY:='5';
               BasePointX:='3';
               BLS24TwistGeneratorX.a.SetFormString('0x5005AD76005C8C92252A3E5320E8CDD88C23900118AAC24EA9DE3020B776C89768182380DDEAD1A259E7723F5D4154FD21A287BE3BE308AB3D9+u*0x431E69C87FF32169AF44356E58175A43B90422F42312250765B6F370A5E7D7D5E17EDECDA87156A81E0B5C2624B3D94410A73A9AB7BCC7F7277');
               BLS24TwistGeneratorX.b.SetFormString('0x8F88538FD214F97EF064012A1A7DDF144AC60BDBBB8B20B5B860AFBED3A889D6EA890CBC6D822FC62EDC4DB924A0AFC24516B49A38CCD034D2+u*0x2C7215B60C3F3D5B1DE11FF681444713B36165186AFA63EFCB9CC2EFEF4C86C954AE6DADB5EC3A2B20715EA44EA7E3031C78F7CDE4D783E91D');
               BLS24TwistGeneratorY.a.SetFormString('0x1C9048A9AF80C00A760219E497D9BE9D17290FC0A7860889BA5E31C7DDD75C8C2259594EADD854095BF6711999F1C1326DD26B4AA61466E707B+u*0x1664F37221CBF9D6C3C5FBBDF09FA299D83257027273EB88272D2A0E73D80C934930E6FE1C603FE4A3BFA7C2769FA92478BA80EF2D8DD87B9EF');
               BLS24TwistGeneratorY.b.SetFormString('0x35D0DF857772F4C7C1B7E51B828679314B22B2A19CBA1D1DE813CB1CFC245145AF18843C18FBD727AE2D7985A4D52B6E675553AD1C11EC2E5ED+u*0xB0C1FA8080D0C0BA67C43AA90EA2FE942DDDF7C269F70BE550A1F8643C23C21547B76BFEDB7E6E4930D83FC4C1D81D2A4240A2C2CA2E2EF8F1');
        end;
scBLS24at192_3:begin
               SecurityLevel:='192bit';
               Family:=cfBLS24;
               u:='0xFFFFC4000000';
               P:='0x55548D56284426D648BC0FF673D986B7C76BA5306446D28ECC59E55F42957C3912A3EE719C7BB39A61C11E3CFFBE150055552D555A05AAAA96AAAAAB';
               new(p.InverseFordivision);
               p.Limit:=2*p.Data.i16[-2];
               p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
               p.id:='Y';
               r:='0xFFFE200189BF476E3612BBDC7C10298561C9BAC29B80FFFF0000EFFFABA00D2EFF3A3F00000000000000000000000001';
               Tr:='0xFFFFC4000001';
               Lp:='0xFFFFC4000000';
               n:='0x55548D56284426D648BC0FF673D986B7C76BA5306446D28ECC59E55F42957C3912A3EE719C7BB39A61C11E3CFFBE150055552D555A04AAAAD2AAAAAB';
               htw:='0x328FE6BCB0BC13F3164C95F0089A7AD6882671535D1BFE0031A6C8CFD60E681E474C495096698F3FA26F09BF6F705699C35F3B7037A32D69697C47'
               +'66A13E03743719076CA3C5E7A78919AAA88DC15EA3D7533CBB9C935129A96EA62BB9CEF8121827B87B81EC558D1741D7D53F8C674986EEF5E3516157BD24'
               +'03B71BC864872E8B6EB058B22BC20174BC9F331EF92695C83BD8E6F4411045D78409AB3E05EA6CA1BF3DD7ADD136A40000F546868E6C1BF7B6088F27A27E6AFABF35B89E06524';
               rtw:=r;
               A:=0;
               B:=5;
               H:='0xFFFFC3FFFFFF';
               FieldParam:=InitFieldParams(-1,p);
               New(FieldParam.MontgomeryData);
               InitMontgomeryStruct(P,FieldParam.MontgomeryData);
               New(TowerParam);
               TowerParam^.FieldParam:=FieldParam;
               TowerParam^.Sigma.SetFieldParams(TowerParam^.FieldParam);
               TowerParam.Sigma.SetFormString('3+u*1');
               TowerParam.pmod8:=p mod 8;
               TwoinFp2.SetFieldParams(TowerParam^.FieldParam);
               TwoinFp2.SetFormString('2+u*0');       // for squaring over FP4
               TowerParam^.inv_2.SetFieldParams(TowerParam^.FieldParam);
               TowerParam^.inv_2:=TwoinFP2.Inverse;
               New(TowerParam2);
               TowerParam2.Tower8paprms:=TowerParam;
               TowerParam2.pmod8:=TowerParam.pmod8;
               TowerParam2.Sigma:=Towerparam.Sigma;
               TowerParam2.FieldParam:=FieldParam;
               TowerParam2.Gamma.SetTowerParams(TowerParam);
               TowerParam2.Gamma.SetFromStrings('0+u*0','1+u*0');
               BtwFp4.SetTowerParams (TowerParam);
               TwistMode:=twDType;
               BtwFp4.a.SetFormString('0+u*0');
               BtwFp4.b.SetFormString('0x2AAA46AB1422136B245E07FB39ECC35BE3B5D29832236947662CF2AFA14ABE1C8951F738CE3DD9CD30E08F1E7FDF0A802AAA96AAAD02D5554B555557+u*0x2AAA46AB1422136B245E07FB39ECC35BE3B5D29832236947662CF2AFA14ABE1C8951F738CE3DD9CD30E08F1E7FDF0A802AAA96AAAD02D5554B555555');
               LoopBin^:=Lp.ToIntArray;
               LoopNaf^:=Lp.ToNafArray;
               LoopTateBin^:=n.ToIntArray;
               LoopTateNaf^:=n.ToNafArray;
               BasePointY:='5';
               BasePointX:='-1';
               BLS24TwistGeneratorX.a.SetFormString('0x121C0CF89301A7319E020D27C8006365F07A3F7313051970E6D2534EA4D302CEE00AC80B96734B8432BABC5FB5CBBA13DF1264D4BD8EB619125F2D9B+u*0x520A509C83836A5AC9BFB37AD53BD33DD66D9E9BA15257CEFCF06CA33594794E6C4C4D4F0AC9AFC9FF1F2A6D8C322D5264B0EEB6C06E2753EF633CBA');
               BLS24TwistGeneratorX.b.SetFormString('0x2050D4E7B2931F8F72A1A3E91F48B654B21137A3D805255C63987DBD2405A58949B5AE61EB9E7901CEC61E790D69A5055120AEA763EDA5FE5DBF1D48+u*0x5009A11B25FFADB3D7CB41B1B5FFBDD641A90BF14BEE58168853FF9A6E0821D7A8352CB8BA5059897B2709753C5E1A5A69FE4BDF1E41986A9008B45B');
               BLS24TwistGeneratorY.a.SetFormString('0x3E97AB88D7D8127D2808C08C987A2EC6979D323DCC540C1A1230A9824B2B733A41D0534D284D1E95D875A5744A9EA7F7BFECC02675D6174D0D09D09C+u*0x422C155C2521DAC82BF07E822AB5B42D151C13FE9DFD86B112736EA42AE1A96C38F2CC28FC32FC0AEEE0343F3E569EF0C7B24B8D1FC5E62D7DA5A4DC');
               BLS24TwistGeneratorY.b.SetFormString('0x29F000AF16684EB37743BA62147EA5BA541C9ADE36F73FC365BE0113705C1D9E63363613A78438C1BC278EF531E98220D155539476A22CEA82931478+u*0x815986A77C2B0C2C6A82BEEF5DDEA4F64BA5E040F249D9774ABC214E416A3642B2C0E468B4B45F8EEF581C988F592A9E12C8547F713BC7801886CC');
        end;
scBLS24at320:begin
               SecurityLevel:='320bit';
               Family:=cfBLS24;
               u:='0x7FFFFFFFFFFFFFFFFC000024000';
               P:='0x15555555555555554EAAAAE6AAA55555564555447557214547C155ACC6AA67AEAE44B685FDF61AE8E4EA557B1B0B0E67F26EBDF13047BD5A7A0090F278562213E93D42C96616'
               +'775C405B3A2C30ABE7910DA7C28B360E5F7F943A57F0B3664F70BFFCE90D7A8B13A0452AF3FAE8E24AC7D1F4768AD3BA5A3F6B23EEED2177749FFFEC5AAB6AAB';
               new(p.InverseFordivision);
               p.Limit:=2*p.Data.i16[-2];
               p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
               p.id:='Y';
               r:='0xFFFFFFFFFFFFFFFFC00002400000000006FFFF82000236FFFF90000BCFFF95B001434FFF6280084E3FCE0E80750C5FA768031DC7F1D4541B4E0BE71440AA6F8D89DF63EFE5'
               +'6C66B514A4B8BF04F1507D7EC00B095C8E584AA43AB4D73F1A000B63FFE65F00000000000001';
               Tr:='0x7FFFFFFFFFFFFFFFFC000024001';
               Lp:='0x7FFFFFFFFFFFFFFFFC000024000';
               n:='0x15555555555555554EAAAAE6AAA55555564555447557214547C155ACC6AA67AEAE44B685FDF61AE8E4EA557B1B0B0E67F26EBDF13047BD5A7A0090F278562213E93D42C9661677'
               +'5C405B3A2C30ABE7910DA7C28B360E5F7F943A57F0B3664F70BFFCE90D7A8B13A0452AF3FAE8E24AC7D1F4768AD3BA5A3F6323EEED217774A0002C5AA92AAB';
               htw:='0x329161F9ADD3C0CA12F68684BD6E9E066AAAA8F1C7552E9C45C733600C861DFB62B0ED3B33BB69FF63383687D29DBCFA94370F5A62FA09FAAECCC0B72A95804485B74CB2920'
               +'E7E45011CBCD9DBA489E21AAC91101F31B36D7A728A7720D77493224FBB53273EE8361B7806A2871A36F1B46D02E6C9F28378D75552870F07B449610BEFA3FC74DC4F731FE2410698'
               +'BF6B86F16D0426F901F110D355A0F80802137FEE6D04FE333CB8DD35273C8B83B4BE843919AB9E1009D791A1A067581D02EBF73D17337A8AD9E40C106776D495EA0AA4021C1FC0D6E'
               +'3301D93CB0CE68FFCD445E36BBE72EF6F820AFF87373DC4A755556DB295E9A9D2502085997E66442D28D424A12612DC2BC53635734798DA7D9C67157CB1C597BB1ED9743EC3108ADC'
               +'987B1BD63232C794AF4EF25CB43CF66C868D4F202725DC5EAE61BB7880F90AE73F4B0D32BA70DCA03D46CBE949E63D8942C1DB2A257FC4C15F7572E76B7AA64F16A94E90FFF4FD5E2'
               +'2A42DF7D1D02EAD27A37AC73A7AEF00CC7D98EBD92CF3C4E05DF3F4F62796575A29476F25ED2C15690CDA0B454698C3CC3F5F241A916327B7303C3E58EE06BB7B81A6524';
               rtw:=r;
               A:=0;
               B:=-2;
               H:='0x7FFFFFFFFFFFFFFFFC000023FFF';
               FieldParam:=InitFieldParams(-1,p);
               New(FieldParam.MontgomeryData);
               InitMontgomeryStruct(P,FieldParam.MontgomeryData);
               New(TowerParam);
               TowerParam^.FieldParam:=FieldParam;
               TowerParam^.Sigma.SetFieldParams(TowerParam^.FieldParam);
               TowerParam.Sigma.SetFormString('1+u*1');
               TowerParam.pmod8:=p mod 8;
               TwoinFp2.SetFieldParams(TowerParam^.FieldParam);
               TwoinFp2.SetFormString('2+u*0');       // for squaring over FP4
               TowerParam^.inv_2.SetFieldParams(TowerParam^.FieldParam);
               TowerParam^.inv_2:=TwoinFP2.Inverse;
               New(TowerParam2);
               TowerParam2.Tower8paprms:=TowerParam;
               TowerParam2.pmod8:=TowerParam.pmod8;
               TowerParam2.Sigma:=Towerparam.Sigma;
               TowerParam2.FieldParam:=FieldParam;
               TowerParam2.Gamma.SetTowerParams(TowerParam);
               TowerParam2.Gamma.SetFromStrings('0+u*0','1+u*0');
               BtwFp4.SetTowerParams (TowerParam);
               TwistMode:=twDType;
               BtwFp4.a.SetFormString('0+u*0');
               BtwFp4.b.SetFormString('0x15555555555555554EAAAAE6AAA55555564555447557214547C155ACC6AA67AEAE44B685FDF61AE8E4EA557B1B0B0E67F26EBDF13047BD5A7A0090'+'F278562213E93D42C96616775C405B3A2C30ABE7910DA7C28B360E5F7F943A57F0B3664F70BFFCE90D7A8B13A0452AF3FAE8E24AC7D1F4768AD3BA5A3F6B23EEED2177749FFFEC5AAB6AAA+u*1');
               LoopBin^:=Lp.ToIntArray;
               LoopNaf^:=Lp.ToNafArray;
               LoopTateBin^:=n.ToIntArray;
               LoopTateNaf^:=n.ToNafArray;
               BasePointY:='5';
               BasePointX:='3';
               BLS24TwistGeneratorX.a.SetFormString('0x1521864DD1D45FD5242F9E9E8E708F5D856419F4C60BF1CA88DF7BCC99BC60B4B4EE94E775427B1AED826A463A83EEFD88433451BACCCEDD988AC3A7E'+'1371AC03B535C1418AA26B870A701DDDE2017F6A1441A6B4D7D27A449461F8149FA6E89E8C4BD36E75BE0BFBD3E04DA2A82483C2ABBF04DBFFE7023FC7F149F74DF147EA15D23E5352B'
                                                  +'+u*0x1DD1E34F7253E3888FA220F57AECAF712B73804E71BA440B3FB41FAF8BBEE4C067BA8957356AC7BC2CE34D97BC64484D2F6FD2D6F8A920835E035324'+'D7579CD54E78CCB6376D9A73C97FD9B8DDE4F21FF36BE44BFAF9503416D7A28F6C5C63ED94DEFA92C49DE4FA3DF086CFBDA63A3E175018FB6241A5C56C2DABC935EF8E261FFA5283705');
               BLS24TwistGeneratorX.b.SetFormString('0x3A56657B056B038481C59724902529B0647EC01C645B421CF32A25A795E998C15995C3B123608B59DA7A921065B72E523AC3BD315A43BE96E07E109733'+'032AFA71DB043F5BE23A31D2AB5FC304B12F78E60DB236189BB3420BB36ABC16ED2BA4E119B42A058022F426EA9B419F57D68CC4411BECD7A024EF8178347181881EF720ED2380E87'
                                                    +'+u*0xC33A8226F22B9C8FBB50030653B010A726E0720B4D713F038D7DDF1BBBDC2D314B9062636612CB6F2AD9972AD7DD9C9F074F477D58C6AB617B44E4'+'DEDA0939F8BC19D1E764D1391631DA57F3467D0D3870991D4FDC41A7524C049AC1751C89BFC28811389A3EDE62C34C24A441E70626674F5C672D05460F9FCCC776FC28797B2534FD45140');
               BLS24TwistGeneratorY.a.SetFormString('0x8883A7FC4C4ED86029F247A172227E292F1ED8565362D45AAB6C9B62001CA13352FE6E55AC6496BE597EB2545D3F60ADD95A09F89F0F2CFD8C4277366A'+'57434E81056E53C1D6EE6422FAD2F71E15816BF8D108E58A2DC0DA04CBB7F57E15F0610DC0A7B1A886DD17F0D2E6BE22ACAE0EC3E58B7438FE1FCFEC62E8414DCDF2493712252B92'
                                                    +'+u*0xD5A28C5229AA9BF953443DE8538FD7CF5A041E7DCBBC6F18B30C5E0A107C133B967A5DE4A11B5697F0361A1D8542CA164ED45C3FEDA58CFF29D4FA'+'6F0CF936B05B87DC3B37E7EBA0CBE54AE7232074F96FD640E01CC7A6BFD259A6000B5AEE35FEFA6A93F1D2822B8203FBC1BB2E87D74C7A7F423B6A25B6D057F84A90EF397AACFFF3A8D64');
               BLS24TwistGeneratorY.b.SetFormString('0xE04B86E8A1979CD9175D59020AAA770CE79C8A7DB372054406FE1F0E979287A3B7EA918E57BE20340C333E06B410910FC16C2BF444B6D711681A0196DA'+'5CE4D9770985EE7CB8315840182FF6544BDAB79746D54A34F1DCC13E2EE83BD06966CBBD5A4AF84534D8E85315CEB2AC2E34FAAB0FB5DF180D3BF2E1F0EE68148A9BB4105E92049FE'
                                                    +'+u*0xEFB9EC45DA46B8B7A441B18BCB3F36BDC471DA52CD24C977286ACAEA52D9F356DB315C0E183629E0826E9E5C5C01EAA49150943EEE8FF0FE54ECCF'+'3A8FE934100057841B0C831EC644787E0BB90B2DEF6DE414210A5C61D2F56C315F5BEDA450C80575C4803148F71BA56DBE5D48617F3CE4CB00ED6F5BCE2198CB6F76005A31D4345EDE2DF');
        end;
end;
TowerParam2^.FrobeniusP_Const_Fp4[0].SetTowerParams(TowerParam);
TowerParam2^.FrobeniusP_Const_Fp4[0]:=TowerParam2^.Gamma.Pow((P-1)/6);
for i:=1 to 4 do begin
                 TowerParam2^.FrobeniusP_Const_Fp4[i].SetTowerParams(TowerParam);
                 TowerParam2^.FrobeniusP_Const_Fp4[i]:=TowerParam2^.FrobeniusP_Const_Fp4[i-1]*TowerParam2^.FrobeniusP_Const_Fp4[0];
                 end;
FrobeniusMapConstX_Fp4.SetTowerParams(TowerParam);
case TwistMode of
twDType: begin
        FrobeniusMapConstX_Fp4:=TowerParam2^.Gamma.Pow((p-1)/3);
        FrobeniusMapConstY_Fp4:=TowerParam2^.Gamma.Pow((p-1)/2);
        end;
twMType:begin
        FrobeniusMapConstX_Fp4:=(TowerParam2^.Gamma.Pow((p-1)/3)).Inverse;
        FrobeniusMapConstY_Fp4:=(TowerParam2^.Gamma.Pow((p-1)/2)).Inverse;
        end;
end;
TowerParam.FrobeniusPi3xSigmaSqr:=(TowerParam2^.Gamma.Pow((p-3)).a);
if TowerParam2.pmod8 mod 4=3 then TowerParam.FrobeniusPi3xSigmaSqr:=TowerParam.FrobeniusPi3xSigmaSqr*TowerParam2.Sigma;
TowerParam2.FrobeniusPi3xSigmaSqr:=TowerParam.FrobeniusPi3xSigmaSqr;
end;
end;



procedure ComputeBLS24Parametres(Params: BLS24CurvesParamsDefinition; var Result:PtrCurveParams);
var localP,LocalR,LocalLp,LocalTr:LInt;
    e,u4,u8,u_12,tmp,tmp1:LInt;
    i:integer;
    _P,_F,m2:HLInt;
    str:String;
    twtmp,twtmp1:Fp4Int;
    tau:array[0..4] of HLint;
    _tau:HLint;
    TwoinFP2:Fp2Int;
    Genrator,Gtmp:Fp4Point;
    BLSTwistBasePointX,BLSTwistBasePointY:Fp4Int;
    FrobeniusPi,FrobeniusPi2:Fp2Int;
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
                Family:=cfBLS24;
                if (Pos('0x',Params.u)<>0)or(Pos('$',Params.u)<>0) then u:=_hex_To_LInt(Params.u)
                else u:=_str_To_LInt(Params.u);
                u4:=(u.Sqr).Sqr;
                u8:=u4.Sqr;
                u_12:=(u-1).sqr;
                r:=u8-u4+1;
                n:=(u_12*r)/3;
                LocalP:=n+u;
                if not IsLIntPrime(LocalP) then raise Exception.Create('Paramétre t invalide pour la construction de la courbe BLS...')
                else begin
                     p:=Localp;
                     p.Limit:=2*p.Data.i16[-2];
                     new(p.InverseFordivision);
                     p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
                     p.id:='Y';
                     Tr:=(u+1);
                     Lp:=u.Absolute;
                     end;
                A:=_Str_To_LInt(Params.A);
                B:=_Str_To_LInt(Params.B);
                New(FieldParam);
                FieldParam.Beta:=Params.Beta;
                FieldParam^.p:=p;
                FieldParam^.p_1_div_2:=(p-1) shr 1;
                FieldParam^.inv_2:=LInt(2).InversModulo(p);
                New(FieldParam.MontgomeryData);
                InitMontgomeryStruct(P,FieldParam.MontgomeryData);
                e:=params.beta;
                if e.IsASqrMod(P,FieldParam.MontgomeryData) then Exception.Create(inttostr(Params.beta)+' is a quadratic Residue modulo P......');

                //computing number of points on E(FP4)with cofactor
                Tau[0]:=2;
                Tau[1]:=tr;
                for i := 1 to 3 do Tau[i+1]:=Tr*Tau[i]-p*Tau[i-1];
                _P:=p.Sqr*p.Sqr;
                _Tau:=Tau[4];
                _F:=(4*_P-_Tau.Sqr) /3;
                _Sqrt_HLint(_F,_F);
                m2:=_P+1-(3*_F+_Tau)/2;
                _VHCopy_LInt(m2,tmp1);
                Htw:=tmp1/r;
                Rtw:=r;         ///  The true order is m2=r*htw , so we can choose a subgourp of order R (the cofactor is htw) (https://eprint.iacr.org/2005/133.pdf)
                H:=(u-1);  // instead of n/r=(u-1)^2/3 , Faster !
                New(TowerParam);
                TowerParam^.FieldParam:=FieldParam;
                TowerParam^.Sigma.SetFieldParams(TowerParam^.FieldParam);
                TowerParam.Sigma.SetFormString(Params.Sigma);
                     ///   should test if beta and sigma are valide parametres
                TwoinFP2.SetFieldParams(FieldParam);
                _Pow_FP2(TowerParam.Sigma,(Result.P.Sqr-1)/2,TwoinFP2);
                if (TwoinFP2.a=1)and(TwoinFP2.b=0) then raise Exception.Create('Paramétre t invalide pour la construction de la courbe BLS...Sigma is a square');
                _Pow_FP2(Result.TowerParam.Sigma,(Result.P.Sqr-1)/3,TwoinFP2);
                if (TwoinFP2.a=1)and(TwoinFP2.b=0) then raise Exception.Create('Paramétre t invalide pour la construction de la courbe BLS...Sigma is a Cube');
                    ///
                TowerParam.pmod8:=p mod 8;
                TwoinFp2.SetFieldParams(TowerParam^.FieldParam);
                TwoinFp2.SetFormString('2+u*0');       // for squaring over FP4
                TowerParam^.inv_2.SetFieldParams(TowerParam^.FieldParam);
                TowerParam^.inv_2:=TwoinFP2.Inverse;
                New(TowerParam2);
                TowerParam2.Tower8paprms:=TowerParam;
                TowerParam2.pmod8:=TowerParam.pmod8;
                TowerParam2.Sigma:=Towerparam.Sigma;
                TowerParam2.FieldParam:=FieldParam;
                TowerParam2.Gamma.SetTowerParams(TowerParam);
                TowerParam2.Gamma.SetFromStrings(Copy(Params.Gamma,2,pos(',',Params.Gamma)-2),Copy(Params.Gamma,pos(',',Params.Gamma)+1,pos(')',Params.Gamma)-pos(',',Params.Gamma)-1));
                BtwFp4.SetTowerParams (TowerParam);
                TwistMode:=Params.TwistMode;
                if TwistMode=twMType then  _mul_fp_fp4(B,TowerParam2.Gamma,BtwFp4) //raise Exception.Create('Invalide parameter   for the generator of the curve, this code supports only D-Type Twistes ....');
                else _mul_fp_fp4(B,TowerParam2.Gamma.Inverse,BtwFp4);
                LoopBin^:=Lp.ToIntArray;
                LoopNaf^:=Lp.ToNafArray;
                LoopTateBin^:=n.ToIntArray;
                LoopTateNaf^:=n.ToNafArray;
                if (Pos('0x',Params.GenratorX)<>0)or(Pos('$',Params.GenratorX)<>0) then BasePointX:=_hex_To_LInt(Params.GenratorX)
                else BasePointX:=_Str_To_LInt(Params.GenratorX);
                ///   Constructing the generator for G1
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
                TwsistGeneratorSeed:=Params.TwistGeneratorSeed;
                BLSTwistBasePointX.SetTowerParams(TowerParam);
                BLSTwistBasePointY.SetTowerParams(TowerParam);
                                ///   Constructing the generator for G2 is done when SetToDefaultGenerator is called
                twTmp.SetTowerParams(TowerParam);
                twTmp1.SetTowerParams(TowerParam);
                BLSTwistBasePointX.a.a.SetToRandom(P,TwsistGeneratorSeed);
                BLSTwistBasePointX.a.b.SetToRandom(P,TwsistGeneratorSeed);
                BLSTwistBasePointX.b.a:=0;
                BLSTwistBasePointX.b.b:=0;
                repeat
                  _Inc_LInt(BLSTwistBasePointX.a.a,1);
                  _Sqr_FP4(BLSTwistBasePointX,twtmp);
                  _Mul_FP4(twtmp,BLSTwistBasePointX,twtmp1);
                  _Add_FP4(twtmp1,BtwFp4,twtmp);
                  until twtmp.IsASquare;
                _Sqrt_FP4(twtmp,BLSTwistBasePointY);

                TowerParam2^.FrobeniusP_Const_Fp4[0].SetTowerParams(TowerParam);
                TowerParam2^.FrobeniusP_Const_Fp4[0]:=TowerParam2^.Gamma.Pow((P-1)/6);
                for i:=1 to 4 do begin
                                 TowerParam2^.FrobeniusP_Const_Fp4[i].SetTowerParams(TowerParam);
                                 TowerParam2^.FrobeniusP_Const_Fp4[i]:=TowerParam2^.FrobeniusP_Const_Fp4[i-1]*TowerParam2^.FrobeniusP_Const_Fp4[0];
                                 end;

                FrobeniusMapConstX_Fp4.SetTowerParams(TowerParam);
                case TwistMode of
                    twDType: begin
                             FrobeniusMapConstX_Fp4:=TowerParam2^.Gamma.Pow((p-1)/3);
                             FrobeniusMapConstY_Fp4:=TowerParam2^.Gamma.Pow((p-1)/2);
                             end;
                twMType:begin
                        FrobeniusMapConstX_Fp4:=(TowerParam2^.Gamma.Pow((p-1)/3)).Inverse;
                        FrobeniusMapConstY_Fp4:=(TowerParam2^.Gamma.Pow((p-1)/2)).Inverse;
                        end;
                end;
                TowerParam.FrobeniusPi3xSigmaSqr:=(TowerParam2^.Gamma.Pow((p-3)).a);
                if TowerParam2.pmod8 mod 4=3 then TowerParam.FrobeniusPi3xSigmaSqr:=TowerParam.FrobeniusPi3xSigmaSqr*TowerParam2.Sigma;
                TowerParam2.FrobeniusPi3xSigmaSqr:=TowerParam.FrobeniusPi3xSigmaSqr;
                Genrator.SetCurveParams(Result);
                Gtmp.SetCurveParams(Result);
                Gtmp.X:=BLSTwistBasePointX;
                Gtmp.Y:=BLSTwistBasePointY;
                Gtmp.Z.a.SetFormString('1');
                Gtmp.Z.b.SetFormString('0');
                Gtmp.Infinity:=false;
                _Mul_Fp4_FpPoint(htw,Gtmp,Genrator);
                BLS24TwistGeneratorX:=Genrator.X;
                BLS24TwistGeneratorY:=Genrator.Y;
                end;
end;

function BLS24CurvesParamsDefinitionToSrtingList(params:BLS24CurvesParamsDefinition):TstringList;
begin
Result:=TStringList.Create;
Result.Add('Parameter u :'+params.u);
if Params.Beta<0 then Result.Add('Irrudictible Polynomial  for Fp2: X^2-'+inttostr(Abs(params.Beta)))
else Result.Add('Irrudictible Polynomial  for Fp2: X^2+'+inttostr(Abs(params.Beta)));
Result.Add('Irrudictible Polynomial  for Fp6: X^3+'+params.Sigma);
if params.B[1]<>'-' then Result.Add('Curve : y^2=x^3+'+params.B) else Result.Add('Curve : y^2=x^3'+params.B);
if params.TwistMode=twDType  then Result.Add('Twiste Mode: D-Type')
else Result.Add('Twiste Mode: M-Type');
Result.Add('Hamming Weight of u (Binary)'+inttostr(HammingWeight(Lint(params.u),false)));
Result.Add('Hamming Weight of u (NAF)'+inttostr(HammingWeight(Lint(params.u),True)));
end;



{ TBLS24Curve }


procedure TBLS24CurvePairing.GenerateParamsTreeView(Tree:TTreeView);
var s:string;
    r:real;
    G:FpPoint;Gtw:Fp4Point;
    tmp:TTreeNode;
begin
Tree.Items.Clear;
tmp:=Tree.Items.Add(nil,'Curve');
Tree.Items.AddChild(tmp,'Family : BLS24');
if CurveParams.B<0 then  s:='y^2=x^3-'+CurveParams.B.Absolute.ToDecimalString
else s:='y^2=x^3+'+CurveParams.B.ToDecimalString;
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
Tree.Items.AddChild(tmp,'Equation of the Twist : '+'y^2=x^3+'+Curveparams.BtwFp4.toHexString);
Tree.Items.AddChild(tmp,'Degree of the Twist : Sextic');
if Curveparams.TwistMode=twMType then S:='M-Type' else S:='D-Type';
Tree.Items.AddChild(tmp,'Type of the Twist : '+S);
tmp:=Tree.Items.Add(nil,'Field G2 (Twist) :E''(Fp4)[R]');
Tree.Items.AddChild(tmp,'Order of G2 R = '+CurveParams.Rtw.ToHexString);
Tree.Items.AddChild(tmp,'Size of G2 : '+inttostr(CurveParams.Rtw.BitLength)+' bit');
Tree.Items.AddChild(tmp,'Cofactor of G2 Htw = '+CurveParams.Htw.ToHexString);
tmp:=Tree.Items.Add(nil,'Targted Field (GT):');
Tree.Items.AddChild(tmp,'Extension Field : Fp24');
Tree.Items.AddChild(tmp,'Size of GT : '+inttostr(24*CurveParams.P.BitLength)+' bit');
tmp:=Tree.Items.Add(nil,'Security Level');
Tree.Items.AddChild(tmp,'Security in G1 : '+inttostr(CurveParams.R.BitLength div 2)+' bit').ImageIndex:=2;
Tree.Items.AddChild(tmp,'Security in G2 : '+inttostr(CurveParams.R.BitLength div 2)+' bit').ImageIndex:=2;
r:=1.92* exp((1/3)*ln(24*CurveParams.P.BitLength))*sqr(exp((1/3)*(ln(ln(24*CurveParams.P.BitLength)/ln(2)))));
Tree.Items.AddChild(tmp,'Security in GT (GNFS) : '+inttostr(Round(r))+' bit').ImageIndex:=2;
tmp:=Tree.Items.Add(nil,'Implemented Pairings');
Tree.Items.AddChild(tmp,'Otp-Ate');
tmp:=Tree.Items.Add(nil,'Tower Construction');
if Curveparams.FieldParam.Beta<0 then Tree.Items.AddChild(tmp,'Fp2<u>=ExstensionField<u,|u^2-'+inttostr(Abs(Curveparams.FieldParam.Beta))+'>')
else Tree.Items.AddChild(tmp,'Fp2<u>=ExstensionField<u,|u^2+'+inttostr(Curveparams.FieldParam.Beta)+'>');
Tree.Items.AddChild(tmp,'Fp4<v>=ExstensionField<v,|v^2-('+Curveparams.TowerParam.Sigma.toHexString+')>');
Tree.Items.AddChild(tmp,'Fp8<w>=ExstensionField<w,|w^2-v>');
Tree.Items.AddChild(tmp,'Fp24<z>=ExstensionField<z,|z^3-w>');
tmp:=Tree.Items.Add(nil,' = '+floattostrf(CurveParams.P.BitLength/CurveParams.R.BitLength,ffGeneral, 2, 4 ));
tmp:=Tree.Items.Add(nil,'Default G1 Generator');
G.SetCurveParams(CurveParams);
G.SetToDefaultGenerator;
Gtw.SetCurveParams (CurveParams);
Gtw.SetToDefaultGenerator;
Tree.Items.AddChild(tmp,G.toHexString);
tmp:=Tree.Items.Add(nil,'Default G2 Generator');
Tree.Items.AddChild(tmp,Gtw.toHexString);
Tree.FullExpand;
Tree.Selected:=Tree.Items[0];
//Tree.SetFocus;
end;

{*******************************************************************************}
{constructor TBLS24CurvePairing.Create(inParametres:StandardBLS24Curves);
begin
CurveParams:=nil;
SetStandardCurveParametres(inParametres);
CoordinatesSystem:=csProjective;
Fp24PoweringMode:=pmNormal;
SetLoopMode(lpmAuto);
PairingAlgorithm:=blspOptAte;
end;          }

{*******************************************************************************}
constructor TBLS24CurvePairing.Create;
begin
inherited create(AOwner);
CurveParams:=nil;
SetStandardCurveParametres(scBLS24at256_1);
CoordinatesSystem:=csProjective;
Fp24PoweringMode:=pmKarbina;
SetLoopMode(lpmAuto);
PairingAlgorithm:=blspOptAte;
end;

{*******************************************************************************}
destructor TBLS24CurvePairing.Destroy;
begin
Dispose(CurveParams);
inherited;
end;

{*******************************************************************************}
function TBLS24CurvePairing.getA: String;
begin
Result:=_A.ToHexString;
end;

{*******************************************************************************}
function TBLS24CurvePairing.getAtw: String;
begin
Result:=CurveParams^.AtwFp4.toHexString;
end;

{*******************************************************************************}
function TBLS24CurvePairing.getB: String;
begin
Result:=_B.ToHexString;
end;

{*******************************************************************************}
function TBLS24CurvePairing.getBeta:integer;
begin
Result:=CurveParams^.FieldParam.Beta;
end;

{*******************************************************************************}
function TBLS24CurvePairing.getBtw: String;
begin
Result:=CurveParams^.BtwFp4.toHexString;
end;

{*******************************************************************************}
function TBLS24CurvePairing.getCoord: CoordinatesSystem;
begin
Result:=CurveParams^.CoordSys;
end;

{*******************************************************************************}
function TBLS24CurvePairing.GetDefautG1Generator: FpPoint;
begin
Result.SetToDefaultGenerator;
end;

{*******************************************************************************}
function TBLS24CurvePairing.GetDefautG2Generator: Fp4Point;
begin
Result.SetToDefaultGenerator;
end;

{*******************************************************************************}
function TBLS24CurvePairing.getGamma: String;
begin
if  TwistMode=twDType then Result:=CurveParams^.TowerParam2.Gamma.toHexString
else Result:=CurveParams^.TowerParam2.Gamma.Inverse.toHexString
end;

{*******************************************************************************}
function TBLS24CurvePairing.getH: string;
begin
Result:=CurveParams^.H.toHexString;
end;

{*******************************************************************************}
function TBLS24CurvePairing.getHtw: string;
begin
Result:=CurveParams^.Htw.toHexString;
end;

{*******************************************************************************}
function TBLS24CurvePairing.getLp: string;
begin
Result:=CurveParams^.Lp.toHexString;
end;

{*******************************************************************************}
function TBLS24CurvePairing.getN: string;
begin
Result:=CurveParams^.N.toHexString;
end;

{*******************************************************************************}
function TBLS24CurvePairing.getP: string;
begin
Result:=CurveParams^.P.toHexString;
end;

{*******************************************************************************}
function TBLS24CurvePairing.getR: string;
begin
Result:=CurveParams^.R.toHexString;
end;

{*******************************************************************************}
function TBLS24CurvePairing.GetRandomG1Point: FpPoint;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsRandomTorsionPoint;
end;

{*******************************************************************************}
function TBLS24CurvePairing.GetRandomG2Point: Fp4Point;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsRandomTorsionPoint;
end;

{*******************************************************************************}
function TBLS24CurvePairing.getRtw: String;
begin
Result:=CurveParams^.Rtw.toHexString;
end;

{*******************************************************************************}
function TBLS24CurvePairing.getSigma: String;
begin
Result:=CurveParams^.TowerParam.Sigma.toHexString;
end;

{*******************************************************************************}
function TBLS24CurvePairing.getTr: string;
begin
Result:=CurveParams^.Tr.toHexString;
end;

{*******************************************************************************}
function TBLS24CurvePairing.GetTwm: TTwistModel;
begin
Result:=CurveParams.TwistMode;
end;

{*******************************************************************************}
function TBLS24CurvePairing.getu: String;
begin
Result:=_u.ToHexString;
end;

{*******************************************************************************}
function TBLS24CurvePairing.HashToG1Point(id: string): FpPoint;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsTorsionFromString(id);
end;

{*******************************************************************************}
function TBLS24CurvePairing.HashToG2Point(id: String): Fp4Point;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsTorsionFromString(id);
end;

{*******************************************************************************}
procedure TBLS24CurvePairing.RecomputeParametres;
var tmp:BLS24CurvesParamsDefinition;
begin
tmp.u:=_u.toHexString;
tmp.Beta:=beta;
tmp.Sigma:=Sigma;
tmp.Gamma:=Gamma;
tmp.A:=_A.toHexString;
tmp.B:=_B.toHexString;
tmp.SecurityLevel:=SecurityLevel;
ComputeBLS24Parametres(tmp,CurveParams);
end;

{*******************************************************************************}
procedure TBLS24CurvePairing.SetStandardCurveParametres(inParametres:StandardBLS24Curves);
begin
case inParametres of
scBLS24at256_0: begin
         FastLoadParams(scBLS24at256_0,CurveParams);
         _A:=BLS24_256_0_Params.A;
         _B:=_Str_To_LInt(BLS24_256_0_Params.B);
          _u:=_Hex_To_LInt(BLS24_256_0_Params.u);
         end;
scBLS24at256_1:begin
         FastLoadParams(scBLS24at256_1,CurveParams);
         _A:=BLS24_256_1_Params.A;
         _B:=_Str_To_LInt(BLS24_256_1_Params.B);
          _u:=_Hex_To_LInt(BLS24_256_1_Params.u);
        end;
scBLS24at256_2:begin
         FastLoadParams(scBLS24at256_2,CurveParams);
         _A:=BLS24_256_2_Params.A;
         _B:=_Str_To_LInt(BLS24_256_2_Params.B);
          _u:=_Hex_To_LInt(BLS24_256_2_Params.u);
        end;
scBLS24Razat256:begin
         ComputeBLS24Parametres(BLS24_Razvan_256_Params,CurveParams);
         _A:=BLS24_256_2_Params.A;
         _B:=_Str_To_LInt(BLS24_Razvan_256_Params.B);
          _u:=_Hex_To_LInt(BLS24_Razvan_256_Params.u);
        end;

scBLS24at192_1:begin
         FastLoadParams(scBLS24at192_1,CurveParams);
         _A:=BLS24_192_1_Params.A;
         _B:=_Str_To_LInt(BLS24_192_1_Params.B);
          _u:=_Hex_To_LInt(BLS24_192_1_Params.u);
        end;
scBLS24at192_2:begin
         FastLoadParams(scBLS24at192_2,CurveParams);
         _A:= BLS24_192_2_Params.A;
         _B:=_Str_To_LInt(BLS24_192_2_Params.B);
          _u:=_Hex_To_LInt(BLS24_192_2_Params.u);
        end;
scBLS24at192_3:begin
         FastLoadParams(scBLS24at192_3,CurveParams);
         _A:= BLS24_192_3_Params.A;
         _B:=_Str_To_LInt(BLS24_192_3_Params.B);
          _u:=_Hex_To_LInt(BLS24_192_3_Params.u);
        end;
scBLS24at192_Raz:begin
          ComputeBLS24Parametres(BLS24_Razvan_192_Params,CurveParams);
         _A:= BLS24_Razvan_192_Params.A;
         _B:=_Str_To_LInt(BLS24_Razvan_192_Params.B);
          _u:=_Hex_To_LInt(BLS24_Razvan_192_Params.u);
        end;
scBLS24at320:begin
         FastLoadParams(scBLS24at320,CurveParams);
         _A:= BLS24_320_1_Params.A;
         _B:=_Str_To_LInt(BLS24_320_1_Params.B);
          _u:=_Hex_To_LInt(BLS24_320_1_Params.u);
        end;

end;
if HammingWeight(Lp,true)<HammingWeight(lp,false) then PerferedPoweringMode:=lpmNaf
else PerferedPoweringMode:=lpmBinary;
end;

{*******************************************************************************}
procedure TBLS24CurvePairing.SetA(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TBLS24CurvePairing.SetB(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TBLS24CurvePairing.SetBeta(Value: Integer);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TBLS24CurvePairing.setCoord(value: CoordinatesSystem);
begin
if value=csJacobian then begin
                         Messagedlg('Jacobian Coordinates systems has not been implemented for BLS24 Curves, it is usless ...',mtwarning,[mbok],0);
                         exit;
                         end;
CurveParams^.CoordSys:=value;
end;

{*******************************************************************************}
procedure TBLS24CurvePairing.SetCustomCurveParametres(inParametres: BLS24CurvesParamsDefinition);
begin
ComputeBLS24Parametres(inParametres,CurveParams);
end;

{*******************************************************************************}
procedure TBLS24CurvePairing.SetGamma(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TBLS24CurvePairing.SetSigma(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TBLS24CurvePairing.Setu(Value: String);
begin
RecomputeParametres;
end;
{********* Pairing function :Compute the pairing according to the specified algorithm***********}
function TBLS24CurvePairing.Paire(P: FpPoint; Q: Fp4Point): Fp24Int;
begin
if (P.Infinity)or(Q.Infinity) then raise Exception.Create('Can''t compute pairing for infinit points ....');
Result:=OptAtePairing(P,Q);
end;

{******************** Set the loop mode of the pairing :Binary / Negatve Representation*************}
procedure TBLS24CurvePairing.SetLoopMode(Value: LoopPoweringMode);
begin
_LoopMode:=Value;
if _LoopMode=lpmAuto then _LoopMode:=PerferedPoweringMode;

case _LoopMode of
lpmBinary:begin
          case PairingAlgorithm of
          blspOptAte:Loop:=CurveParams.LoopBin;
          end;
          end;
lpmNaf:begin
       case PairingAlgorithm of
       blspOptAte:Loop:=CurveParams.LoopNaf;
       end;
       end;
end;
end;

procedure TBLS24CurvePairing.SetPArams(const Value: StandardBLS24Curves);
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

{******************** Final Exponentiation Step :Compute f^((p^24-1)/r) ***********************}
Procedure TBLS24CurvePairing.FinalPowerOptimalAteBLS1(f:Fp24Int;u:LInt;PoweringMode:FpPoweringMode;var Result:Fp24Int);
var mu :array[0..7]of Fp24Int;
    m1,m2,ft1,ft,mu7i:Fp24Int;
begin
///***** Soft Exponentiation ***///
_Inv_FP24(f,ft);
_Conjugate_FP24(f,ft1);
_Mul_FP24(ft1,ft,ft);
_Pow_FP24_P_i(ft,4,ft1);
_Mul_FP24(ft1,ft,ft);
///***** Hard Exponentiation ***///
_Pow_FP24(ft,u,m1,PoweringMode);
_Pow_FP24(m1,u,m2,PoweringMode);
_Unitary_Sqr_FP24(m1,m1);
_Conjugate_FP24(m1,m1);
_Mul_FP24(m2,m1,mu[7]);
_Mul_FP24(mu[7],ft,mu[7]);
_Pow_FP24(mu[7],u,mu[6],PoweringMode);
_Pow_FP24(mu[6],u,mu[5],PoweringMode);
_Pow_FP24(mu[5],u,mu[4],PoweringMode);
_Conjugate_FP24(mu[7],mu7i);
_Pow_FP24(mu[4],u,mu[3],PoweringMode);
_Mul_FP24(mu[3],mu7i,mu[3]);
_Pow_FP24(mu[3],u,mu[2],PoweringMode);
_Pow_FP24(mu[2],u,mu[1],PoweringMode);
_Pow_FP24(mu[1],u,mu[0],PoweringMode);
_Unitary_Sqr_FP24(ft,ft1);
_Mul_FP24(mu[0],ft1,mu[0]);
_Mul_FP24(mu[0],ft,mu[0]);
_Pow_FP24_P_i(mu[1],1,m2);
_Pow_FP24_P_i(mu[1],1,mu[1]);
_Pow_FP24_P_i(mu[2],2,mu[2]);
_Pow_FP24_P_i(mu[3],3,mu[3]);
_Pow_FP24_P_i(mu[4],4,mu[4]);
_Pow_FP24_P_i(mu[5],5,mu[5]);
_Pow_FP24_P_i(mu[6],6,mu[6]);
_Pow_FP24_P_i(mu[7],7,mu[7]);
_Mul_FP24(mu[0],mu[1],Result);
_Mul_FP24(Result,mu[2],Result);
_Mul_FP24(Result,mu[3],Result);
_Mul_FP24(Result,mu[4],Result);
_Mul_FP24(Result,mu[5],Result);
_Mul_FP24(Result,mu[6],Result);
_Mul_FP24(Result,mu[7],Result);
end;

{******************** Compute the Optimal Ate Pairing ***********************}
Function TBLS24CurvePairing.OptAtePairing(Pt:FpPoint;Qt:Fp4Point):Fp24Int;
var T,mQt:Fp4Point;
    i:integer;
    f,f1,f2:Fp24Int;
begin

f.SetTowerParams(CurveParams.TowerParam2);
f.SetToOne;
T:=Qt;
_Neg_Fp4_Point(Qt,mQt);
T.SetPairingPointCoordinates(Pt.X,Pt.Y);
T.ComputeLigneValue:=true;
for i:=Length(Loop^)-2 Downto 0 do begin
                                   case CoordinatesSystem of
                                      csAffine:_Double_Affine_Fp4_Point(T,T,true);
                                      csProjective:_Double_Projective_Fp4_Point(T,T,True);
                                   end;
                                   _Sqr_FP24(f,f);
                                   _Sparse_Mul_FP24(f,T.LineAtP24^,f,TwistMode);
                                   if  Loop^[i]=1 then begin
                                                        case CoordinatesSystem of
                                                          csAffine:_Add_Affine_Fp4_Point(Qt,T,T,True);
                                                          csProjective:_Add_Projective_Fp4_Point(Qt,T,T,True);
                                                          end;
                                                       _Sparse_Mul_FP24(f,T.LineAtP24^,f,TwistMode);
                                                       end
                                   else if (Loop^[i]=-1) then begin
                                                                case CoordinatesSystem of
                                                                  csAffine:_Add_Affine_Fp4_Point(mQt,T,T,True);
                                                                  csProjective:_Add_Projective_Fp4_Point(mQt,T,T,True);
                                                                  end;
                                                              _Sparse_Mul_FP24(f,T.LineAtP24^,f,TwistMode);
                                                              end;
                                   end;
FinalPowerOptimalAteBLS1(f,Qt.CurveParams.u,Fp24PoweringMode,Result);
end;

procedure Register;
begin
  RegisterComponents('Pairings', [TBLS24CurvePairing]);
end;

end.
