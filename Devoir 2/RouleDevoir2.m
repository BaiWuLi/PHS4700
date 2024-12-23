format long;
clear;
global longueur largeur hauteur hfilet lfilet dir ;
global RayonBalle MasseBalle ;
longueur=274/100;
largeur=152.5/100;
hauteur=76/100;
hfilet=15.25/100;
lfilet=183/100;
for vit=1:4
  if vit == 1
% OK
    rbi=[0.0; 0.5; 1.1];
    vbi=[4; 0; 0.8]; %m/s
    wbi=[0;  -70; 0]; %rad/s
  elseif vit == 2
% OK
    rbi=[0.0; 0.40; 1.14];
    vbi=[10.0; 1.0; 0.2]; %m/s
    wbi=[0; 100; -50]; %rad/s
  elseif vit == 3
% OK
    rbi=[2.74; 0.5; 1.14];
    vbi=[-5.00; 0; 0.2]; %m/s
    wbi=[0;  100; 0]; %rad/s
  elseif vit == 4
% OK
    rbi=[0.0; 0.3; 1.00];
    vbi=[10.0; -2.0; 0.2]; %m/s
    wbi=[0;  10; -100]; %rad/s
  end;
  if rbi(1) < longueur/2-RayonBalle
    dir=1;
  else
    dir=-1;
  end;
%
%  Tracer de la table
%
  clf;
  hold;
  fprintf('\nSimulation %3d\n',vit);
  xlabel('x(m)');
  ylabel('y(m)');
  zlabel('z(m)');
  axis equal;
  TabledePingPong;
  for option=1:3
    if option == 1
      fprintf('Acceleration gravitationnelle seulement \n');
      type='r-';
    elseif option == 2
      fprintf('Acceleration gravitationnelle et force visqueuse \n');
      type='b-';
    elseif option == 3
      fprintf('Acceleration gravitationnelle, force visqueuse et force de Magnus\n');
      type='k-';
    end;
    fprintf('Position initiale de la balle [%12.8f %12.8f %12.8f]  m     \n',rbi(1),rbi(2),rbi(3));
    fprintf('Vitesse initiale de la balle  [%12.8f %12.8f %12.8f]  m/s   \n',vbi(1),vbi(2),vbi(3));
    fprintf('Vitesse angulaire de la balle [%12.8f %12.8f %12.8f]  rad/s \n',wbi(1),wbi(2),wbi(3));
    [coup vbf tt xx yy zz]=Devoir2(option,rbi,vbi,wbi);
    sz=size(tt);
    plot3(xx,yy,zz,type);
    fprintf('\n La simualtion se termine au temps %12.8f s\n',tt(sz));
    if coup == 0
      fprintf('Le coup est reussi \n');
    elseif coup == 1
      fprintf('Le balle touche la table du mauvais cote \n');
    elseif coup == 2
      fprintf('Le balle touche le filet \n');
    elseif coup == 3
      fprintf('Le balle touche le sol\n');
    end
    fprintf('Vitesse finale de la balle     [%12.8f %12.8f %12.8f]  m/s\n',vbf(1),vbf(2),vbf(3));
    fprintf('Position finale de la balle    [%12.8f %12.8f %12.8f]  m \n',xx(sz),yy(sz),zz(sz));
    fprintf('\n\n');
    pause;
  end
  title('Option 1 -> trajectoire en rouge, Option 2 -> trajectoire en bleu et Option 3 -> trajectoire en noir','FontSize',14)
  hold;
end
