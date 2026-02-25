

classdef TwaveLoader
    properties
        data_path_and_ranges = {
          "data/VMedians/Normal.ecg_1.vmedianecg", 498, 734;  % path, begin, end
          "data/VMedians/Alvale12_normal_elec.ecg_1.vmedianecg", 520, 800;
          "data/VMedians/Alvale15_normal_elec.ecg_1.vmedianecg", 545, 800; % P + T wave
          "data/VMedians/Alvale15_normal_elec.ecg_2.vmedianecg", 511, 785;
          "data/VMedians/HCM 1.xml_0.vmedianecg", 550, 930;
          "data/VMedians/HCM 1.xml_1.vmedianecg", 525, 904;
          "data/VMedians/HCM 2.xml_0.vmedianecg", 572, 825;
          "data/VMedians/HCM 2.xml_1.vmedianecg", 506, 775;
          "data/VMedians/HCM 3.xml_0.vmedianecg", 529, 840;
          "data/VMedians/HCM 3.xml_1.vmedianecg", 485, 787;
          "data/VMedians/HCM 4.xml_0.vmedianecg", 537, 837;
          "data/VMedians/HCM 4.xml_1.vmedianecg", 504, 790;
          "data/VMedians/LAD Proximal.ecg_1.vmedianecg", 514, 765;
          "data/VMedians/LBBB male 70 id-4331.ecg_1.vmedianecg", 590, 820;
          "data/VMedians/LBBB male 70 id-8332.ecg_1.vmedianecg", 570, 880; % 15?
          "data/VMedians/LCX 1.ecg_1.vmedianecg", 506, 815;
          "data/VMedians/LCX 2.ecg_1.vmedianecg", 550, 825;
          "data/VMedians/LCX 3.ecg_1.vmedianecg", 555, 811;
          % "data/VMedians/LCX 4.ecg_1.vmedianecg", 711, 860; % 19? Twave neg???
          "data/VMedians/Normal.ecg_1.vmedianecg", 499, 734;
          "data/VMedians/PVC 15f.ecg_1.vmedianecg", 490, 775;
          "data/VMedians/PVC 15f.ecg_2.vmedianecg", 534, 785;
          "data/VMedians/PVC 15f.ecg_3.vmedianecg", 548, 785;
          "data/VMedians/PVC00.ecg_1.vmedianecg", 511, 890;
          % "data/VMedians/PVC00.ecg_2.vmedianecg", 560, 860; % 25? PVC
          "data/VMedians/PVC15a.ecg_1.vmedianecg", 538, 785;
          "data/VMedians/PVC15a.ecg_2.vmedianecg", 557, 830;
          "data/VMedians/PVC34.ecg_1.vmedianecg", 516, 800;
          "data/VMedians/PVC34.ecg_2.vmedianecg", 567, 837;
          "data/VMedians/PVC34.ecg_3.vmedianecg", 555, 820;
          "data/VMedians/PVC41.ecg_2.vmedianecg", 597, 846;
          "data/VMedians/RBBB female 65 id-12552.ecg_1.vmedianecg", 545, 845;
          "data/VMedians/RBBB female 69 id-5451.ecg_1.vmedianecg", 545, 805;
          % "data/VMedians/RCA Distal.ecg_1.vmedianecg", 484, 736; % 34? ST elev
          "data/VMedians/RCA Proximal.ecg_1.vmedianecg", 500, 786;
          "data/VMedians/TYPE 1 N.11 a baseline.zip.ecg_1.vmedianecg", 600, 820;
          "data/VMedians/UTSW_1_ecg1.xml_1.vmedianecg", 500, 695;
          "data/VMedians/UTSW_1_ecg2.xml_1.vmedianecg", 170, 285;
          "data/VMedians/UTSW_1_ecg2.xml_1.vmedianecg", 540, 665;
          "data/VMedians/UTSW_1_ecg3.xml_1.vmedianecg", 100, 270;
          "data/VMedians/UTSW_1_ecg3.xml_1.vmedianecg", 545, 690;
          "data/VMedians/UTSW_1_ecg3.xml_1.vmedianecg", 960, 1100;
          % "data/VMedians/UTSW_1_ecg4.xml_1.vmedianecg", 116, 280; %
          % baseline needs to be fixed for this case
          "data/VMedians/UTSW_1_ecg4.xml_1.vmedianecg", 484, 655;
          % "data/VMedians/UTSW_1_ecg5.xml_1.vmedianecg", 545, 700;
          "data/VMedians/UTSW_1_ecg5.xml_1.vmedianecg", 1020, 1195;
          "data/VMedians/UTSW_1_ecg6.xml_1.vmedianecg", 190, 310;
          "data/VMedians/UTSW_1_ecg6.xml_1.vmedianecg", 550, 650;
          "data/VMedians/UTSW_1_ecg6.xml_1.vmedianecg", 880, 1020;
          "data/VMedians/UTSW_1_ecg7.xml_1.vmedianecg", 530, 670;
          "data/VMedians/UTSW_1_ecg7.xml_1.vmedianecg", 970, 1115;
          "data/VMedians/UTSW_1_ecg8.xml_1.vmedianecg", 550, 710;
          "data/VMedians/UTSW_1_ecg8.xml_1.vmedianecg", 520, 700;
          "data/VMedians/UTSW_1_ecg9.xml_2.vmedianecg", 17, 227;
          "data/VMedians/UTSW_1_ecg9.xml_2.vmedianecg", 490, 690;
          "data/VMedians/UTSW_1_ecg9.xml_2.vmedianecg", 1000, 1180;
          "data/VMedians/UTSW_1_ecg10.xml_1.vmedianecg", 511, 707;
          "data/VMedians/UTSW_1_ecg11.xml_1.vmedianecg", 488, 715;
          "data/VMedians/UTSW_1_ecg12.xml_1.vmedianecg", 510, 725;
          "data/VMedians/UTSW_1_ecg13.xml_1.vmedianecg", 506, 750;
          "data/VMedians/UTSW_1_ecg14.xml_1.vmedianecg", 492, 734;
          "data/VMedians/UTSW_1_ecg15.xml_1.vmedianecg", 530, 730;
          "data/VMedians/UTSW_1_ecg16.xml_1.vmedianecg", 525, 725;
          "data/VMedians/UTSW_1_ecg17.xml_1.vmedianecg", 490, 725;
          "data/VMedians/UTSW_1_ecg18.xml_1.vmedianecg", 490, 745;
          "data/VMedians/UTSW_1_ecg19.xml_1.vmedianecg", 500, 730;
          "data/VMedians/UTSW_1_ecg20.xml_1.vmedianecg", 510, 741;
          "data/VMedians/UTSW_1_ecg21.xml_1.vmedianecg", 500, 730;
          "data/VMedians/UTSW_1_ecg22.xml_1.vmedianecg", 500, 740;
          "data/VMedians/UTSW_1_ecg23.xml_1.vmedianecg", 490, 740;
          "data/VMedians/UTSW_1_ecg24.xml_1.vmedianecg", 500, 740;
          "data/VMedians/UTSW_807_ecg103874104642.xml_1.vmedianecg", 516, 806;
          "data/VMedians/UTSW_807_ecg18781113706.xml_1.vmedianecg", 490, 820;
          % "data/VMedians/UTSW_807_ecg24941101813.xml_1.vmedianecg", 1, 165;
          "data/VMedians/UTSW_807_ecg24941101813.xml_1.vmedianecg", 485, 730;
          % "data/VMedians/UTSW_807_ecg26995184837.xml_1.vmedianecg", 1, 195;
          "data/VMedians/UTSW_807_ecg26995184837.xml_1.vmedianecg", 490, 720;
          "data/VMedians/UTSW_807_ecg27050122825.xml_1.vmedianecg", 90, 290;
          "data/VMedians/UTSW_807_ecg27050122825.xml_1.vmedianecg", 478, 700;
          "data/VMedians/UTSW_807_ecg27050122825.xml_1.vmedianecg", 895, 1100;
          "data/VMedians/UTSW_807_ecg27488085335.xml_1.vmedianecg", 490, 760;
          "data/VMedians/UTSW_807_ecg40657123934.xml_1.vmedianecg", 500, 770;
          "data/VMedians/UTSW_807_ecg51800090458.xml_1.vmedianecg", 500, 780;
          "data/VMedians/UTSW_807_ecg62149102258.xml_1.vmedianecg", 510, 780;
          "data/VMedians/UTSW_807_ecg67323085649.xml_1.vmedianecg", 510, 800;
          "data/VMedians/UTSW_807_ecg72306093634.xml_1.vmedianecg", 510, 790;
          "data/VMedians/UTSW_807_ecg82272094301.xml_1.vmedianecg", 510, 810;
          "data/VMedians/UTSW_807_ecg8706102829.xml_1.vmedianecg", 490, 675;
          "data/VMedians/UTSW_807_ecg88815090125.xml_1.vmedianecg", 500, 790;
          "data/VMedians/UTSW_807_ecg93607101400.xml_1.vmedianecg", 500, 800;
          "data/VMedians/WPW 1.xml_0.vmedianecg", 580, 860;
          "data/VMedians/WPW 1.xml_1.vmedianecg", 500, 790;
          "data/VMedians/WPW 2.xml_0.vmedianecg", 570, 890;
          "data/VMedians/WPW 2.xml_1.vmedianecg", 540, 855;
          "data/VMedians/X006.inf.ecg_1.vmedianecg", 195, 390;
          "data/VMedians/X006.inf.ecg_1.vmedianecg", 540, 735;
          "data/VMedians/X006.inf.ecg_1.vmedianecg", 900, 1060;
          "data/VMedians/X007.inf.ecg_1.vmedianecg", 166, 355;
          "data/VMedians/X007.inf.ecg_1.vmedianecg", 516, 722;
          "data/VMedians/X007.inf.ecg_1.vmedianecg", 883, 1088;
          % "data/VMedians/test_ECGexpert_2025-03-04T11_38_13.ecg_1.vmedianecg", 510, 753;
          "data/VMedians/test_ECGexpert_2025-03-04T11_39_03.ecg_1.vmedianecg", 528, 730;
          "data/VMedians/test_ECGexpert_2025-03-04T11_49_50.ecg_1.vmedianecg", 510, 744;
          "data/VMedians/test_ECGexpert_2025-03-04T12_04_32.ecg_1.vmedianecg", 540, 800;
          "data/VMedians/test_ECGexpert_2025-03-04T12_06_09.ecg_1.vmedianecg", 530, 785;
        };
        data_cnt = -1;
    end
    
    methods
        function obj = TwaveLoader(varargin)
            if(obj.data_cnt < 0)
                obj.data_cnt = length(obj.data_path_and_ranges);
            end
            if nargin > 0
                for k = 1:2:length(varargin)
                    obj.(varargin{k}) = varargin{k+1};
                end
            end
        end
        
        function [twave, tdom] = getVMediansTwave(obj, number)
            ecg = loadmat(obj.data_path_and_ranges{number, 1});
            signal = rms([ecg(4:6,:)/1.5; ecg(7:end,:)]);
            
            twave_start = obj.data_path_and_ranges{number, 2};
            twave_end   = obj.data_path_and_ranges{number, 3};
            twave = signal(1,twave_start:twave_end);

            tdom = 1 - cumsum(twave)/max(cumsum(twave));
        end

        function [name, twave_begin, twave_end] = getVMediansPathAndRanges(obj, number)
            name            = obj.data_path_and_ranges{number, 1};
            twave_begin     = obj.data_path_and_ranges{number, 2};
            twave_end       = obj.data_path_and_ranges{number, 3};
        end

        function signal = getVMediansOriginalSignalsRms(obj, number)
            ecg = loadmat(obj.data_path_and_ranges{number, 1});
            signal = rms([ecg(4:6,:)/1.5; ecg(7:end,:)]);
        end
    end
end

