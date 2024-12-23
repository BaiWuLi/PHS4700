clear
%
% Procedure servant à rouler le devoir 1
%
% Cas 1
AngRotCas1=0.0;
vangulaireCas1=[0;0;0];
forcesCas1=[11e6; 8.75e6; 8.75e6];
posNLCas1=[0;0;0];
[pcmNLCas1 INLCas1 alphaNLCas1]=Devoir1(AngRotCas1,vangulaireCas1,forcesCas1,posNLCas1);
fprintf('\nResultats navette verticale \n');
fprintf('  Centre de masse = ( %10.5f  %10.5f  %10.5f )\n',pcmNLCas1(1),pcmNLCas1(2),pcmNLCas1(3));
fprintf('  Moment inertie  =\n   %14.0f & %14.0f & %14.0f \\\\ \n   %14.0f & %14.0f & %14.0f \\\\ \n   %14.0f & %14.0f & %14.0f \\\\ \n',...
     INLCas1(1,1),INLCas1(1,2),INLCas1(1,3),INLCas1(2,1),INLCas1(2,2),INLCas1(2,3),INLCas1(3,1),INLCas1(3,2),INLCas1(3,3));
fprintf('  acc angulaire   = ( %10.5f  %10.5f  %10.5f )\n',alphaNLCas1(1),alphaNLCas1(2),alphaNLCas1(3));
%
% Cas 2
AngRotCas2=-pi/3;
vangulaireCas2=[-0.54; 0; 00];
forcesCas2=[11e6; 8.75e6; 0];
posNLCas2=[0;-19.6075;50];
[pcmNLCas2 INLCas2 alphaNLCas2]=Devoir1(AngRotCas2,vangulaireCas2,forcesCas2,posNLCas2);
fprintf('\nResultats navette incline\n');
fprintf('  Centre de masse = ( %10.5f  %10.5f  %10.5f )\n',pcmNLCas2(1),pcmNLCas2(2),pcmNLCas2(3));
fprintf('  Moment inertie  =\n   %14.0f & %14.0f & %14.0f \\\\ \n   %14.0f & %14.0f & %14.0f \\\\ \n   %14.0f & %14.0f & %14.0f \\\\ \n',...
     INLCas2(1,1),INLCas2(1,2),INLCas2(1,3),INLCas2(2,1),INLCas2(2,2),INLCas2(2,3),INLCas2(3,1),INLCas2(3,2),INLCas2(3,3));
fprintf('  acc angulaire   = ( %10.5f  %10.5f  %10.5f )\n',alphaNLCas2(1),alphaNLCas2(2),alphaNLCas2(3));
%
% Ancien Cas2
%AngRotCas2=-pi/9;
%vangulaireCas2=[-0.6; 0; 0.1];
%forcesCas2=[11e6; 8.75e6; 0];
%posNLCas2=[0;0;1000];
%[pcmNLCas2 INLCas2 alphaNLCas2]=Devoir1(AngRotCas2,vangulaireCas2,forcesCas2,posNLCas2);
%fprintf('\n  Resultats navette incline \n');
%fprintf('  Centre de masse = ( %10.5f  %10.5f  %10.5f )\n',pcmNLCas2(1),pcmNLCas2(2),pcmNLCas2(3));
%fprintf('  Moment inertie  =\n   %14.0f & %14.0f & %14.0f \\\\ \n   %14.0f & %14.0f & %14.0f \\\\ \n   %14.0f & %14.0f & %14.0f \\\\ \n',...
%     INLCas2(1,1),INLCas2(1,2),INLCas2(1,3),INLCas2(2,1),INLCas2(2,2),INLCas2(2,3),INLCas2(3,1),INLCas2(3,2),INLCas2(3,3));
%fprintf('  acc angulaire   = ( %10.5f  %10.5f  %10.5f )\n',alphaNLCas2(1),alphaNLCas2(2),alphaNLCas2(3));
