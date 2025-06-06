%%
%Energy control and design optimization of a hybrid solar-hydrogen energy storage system using various solar panel technologies
% Journal of Energy Storage Volume 94, 30 July 2024, 112389
%CO2 EMISSIONS OF SILICON PHOTOVOLTAIC MODULES- IMPACT OF MODULE DEISGN AND PRODUCTION LOCATION: 480 kg CO2/kW,
%  17.6 gCO2/kWh
%%

function MODULE_LIB = getPVModuleLibrary()
    MODULE_LIB = struct(...
        'Mono_SI', struct('Embedded_CO2', 17.6,'Degrading_Factor', 0.88, 'panel_efficiency', 0.20,  'temp_coeff', -0.0031,   'ref_temp', 25, 'ncot', 46, 'life', 20, 'IC', 210, 'OMC', 4.2, 'RC', 210), ...
       'Mono_SI_SOLARIA', struct('Embedded_CO2', 17.6,'Degrading_Factor', 0.88, 'ratedpowerkw',0.4,'SurfaceAream2',1.92,'panel_efficiency', 0.2020,  'temp_coeff', -0.0039,   'ref_temp', 25, 'ncot', 45, 'life', 30, 'IC', 145, 'OMC', 145*0.02, 'RC', 145*0.8), ...
...
       'Poly_SI', struct('Embedded_CO2',  17.6,'Degrading_Factor', 0.88, 'panel_efficiency', 0.15,  'temp_coeff', -0.0031, 'ref_temp', 25,'ncot', 46, 'life', 20, 'IC', 120, 'OMC', 2.4, 'RC', 120), ...
        'Thin_Film', struct('Embedded_CO2',  17.6,'Degrading_Factor', 0.88, 'panel_efficiency', 0.07, 'temp_coeff', -0.0031, 'ref_temp', 25,'ncot', 33, 'life', 20, 'IC', 46, 'OMC', 0.92, 'RC', 46) ...
    );
end
%%
%COST of soalr panel
%https://www.irena.org/-/media/Files/IRENA/Agency/Publication/2024/Sep/IRENA_Renewable_power_generation_costs_in_2023.pdf
%Technical Specifications
%Solaria PowerXT 430R-PL https://static1.squarespace.com/static/5fa44724b1ca66220ca21c60/t/627962fbf47df34a3d5c34b2/1652122364846/Datasheet_PowerXT_Resi-400PM.pdf


%https://venturama-solar.de/produkt/ja-solar-jam54s31-lr-420w-full-black-pv-modul/?utm_source=bing&utm_medium=cpc&utm_campaign=DE%20%7C%20PMax%20%7C%20Shop%20%7C%20PV%20Module%20%7C%201700%25&utm_term=2333919589172141&utm_content=Solarmodule%20(PV)