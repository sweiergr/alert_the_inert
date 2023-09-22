%{
    This is just an auxiliary file that allows waf to generate a complete DAG, if you prefer to run estimate_model.m manually and potentially in several chunks to distribute the computational burden.

    This file is not needed if you let waf run estimate_model.m

%}

clc
clear
fprintf('Placebo file to close DAG for WAF is run... Proceeding to counterfactuals and formatting results...\n');