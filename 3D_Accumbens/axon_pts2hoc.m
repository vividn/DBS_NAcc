function [neuron] = axon_pts2hoc(axonpts,somapts,DendStructs)

% axon_pts2hoc()
%
% Matt Johnson (June 2007, revised September 2011)
% Nate Faber (revised February 2014)
%
% Program that reads the XYZ coordinates of a neuron reconstruction from
% Rhino and then creates a series of text files containing .hoc code of the 
% reconstructed neuronal segments.  Inputs are the branched neuron
% structure:
%     neuron - (3 x curvelength x numberofcurves)
%
% if somapts is supplied, the program also writes a create soma line
% and adds 3d points for the soma, and connects it to the initialsegment

h = waitbar(0,'Please wait...');
nNeurons = size(axonpts,3);

% cd ML_hocpoints;

% Parse through axonal tree creating pre-defined segment lengths
figure; hold on;
for iNeuron = 1:nNeurons,

    waitbar(iNeuron/nNeurons,h,['Parsing fiber: ',num2str(iNeuron)]);
    
    neuron(iNeuron).axonseg = axonpointparse(squeeze(axonpts(:,:,iNeuron))*1e3);
    axonseg = neuron(iNeuron).axonseg;
    
    % Calculate numbers for header file
    lastnode = find(axonseg(:,7)==1,1,'last');
    if lastnode == length(axonseg(:,7)),
        lastnode = lastnode - 8;
    end
    plot3(axonseg(1:lastnode,1),axonseg(1:lastnode,2),axonseg(1:lastnode,3),'ko-');
    
    axonnodes = length(find(axonseg(1:lastnode,7)==1));
    paranodes1 = length([find(axonseg(1:lastnode,7)==2);...
        find(axonseg(1:lastnode,7)==8)]);
    paranodes2 = length([find(axonseg(1:lastnode,7)==3);...
        find(axonseg(1:lastnode,7)==7)]);
    axoninter = length([find(axonseg(1:lastnode,7)==4);...
        find(axonseg(1:lastnode,7)==5);...
        find(axonseg(1:lastnode,7)==6)]);
    total = axonnodes + paranodes1 + paranodes2 + axoninter;
    
    % Write the header files
    fid = fopen(['..\axonpts\axon',num2str(iNeuron),'.txt'],'w');
    fprintf(fid,'axonnodes = %d\n',axonnodes);
    fprintf(fid,'paranodes1 = %d\n',paranodes1);
    fprintf(fid,'paranodes2 = %d\n',paranodes2);
    fprintf(fid,'axoninter = %d\n',axoninter);
    fprintf(fid,'totalAxon = %d\n',total);
    fprintf(fid,'nDends = %d\n\n',size(DendStructs,1));

    % Write the create statements
    if nargin >= 2
        fprintf(fid,'create soma\n');
    end
    if nargin >= 3
        nDends = size(DendStructs,1);
        fprintf(fid,'create proxdend[%d], middend[%d], distdend[%d]\n',nDends,nDends*2,nDends*2*2);
    end
    fprintf(fid,'create node[axonnodes], MYSA[paranodes1]\n');
    fprintf(fid,'create FLUT[paranodes2], STIN[axoninter]\n\n');    
    
    % Connection statements
    fprintf(fid,'//connect everything together\n');
    if nargin >= 2
        fprintf(fid,'connect node[0](0), soma(1)\n\n');
    end
    fprintf(fid,[...
        'for i=0, axonnodes-2 {\n',...
        'connect MYSA[2*i](0), node[i](1)\n',...
		'connect FLUT[2*i](0), MYSA[2*i](1)\n',...
		'connect STIN[3*i](0), FLUT[2*i](1)\n',...
		'connect STIN[3*i+1](0), STIN[3*i](1)\n',...
		'connect STIN[3*i+2](0), STIN[3*i+1](1)\n',...
		'connect FLUT[2*i+1](0), STIN[3*i+2](1)\n',...
		'connect MYSA[2*i+1](0), FLUT[2*i+1](1)\n',...
		'connect node[i+1](0), MYSA[2*i+1](1)\n',...
        '}\n\n']);
    if nargin >= 3
        for iDend = 0:nDends-1
            if (iDend+1)/nDends <= 0.5
                somaSide = 1;
            else
                somaSide = 0;
            end
            fprintf(fid,'connect proxdend[%d], soma(%d)\n',iDend,somaSide);
            fprintf(fid,'connect middend[%d](0), proxdend[%d](1)\n',iDend*2,iDend);
            fprintf(fid,'connect middend[%d](0), proxdend[%d](1)\n',iDend*2+1,iDend);
            fprintf(fid,'connect distdend[%d](0), middend[%d](1)\n',iDend*4,iDend*2);
            fprintf(fid,'connect distdend[%d](0), middend[%d](1)\n',iDend*4+1,iDend*2);
            fprintf(fid,'connect distdend[%d](0), middend[%d](1)\n',iDend*4+2,iDend*2+1);
            fprintf(fid,'connect distdend[%d](0), middend[%d](1)\n',iDend*4+3,iDend*2+1);
            fprintf(fid,'\n');
        end
    end
    

    
    % Write 3d points
    fprintf(fid,'//****************\n');
    fprintf(fid,'//COMPARTMENT LIST\n');
    fprintf(fid,'//****************\n\n');
    fprintf(fid,'//beginpoints\n')

    if nargin >= 2
        fprintf(fid,'soma{\n');
        pos1 = somapts(1,:)*1e3;
        pos2 = somapts(2,:)*1e3;
        fprintf(fid,'pt3dadd(%f,%f,%f,0)\n',pos1(1),pos1(2),pos1(3));
        fprintf(fid,'pt3dadd(%f,%f,%f,0)\n',pos2(1),pos2(2),pos2(3));
        fprintf(fid,'}\n\n');
    end
    
    if nargin >= 3
        for iDend = 0:nDends-1
            pDend = DendStructs{iDend+1,nNeurons};
            fprintf(fid,'proxdend[%d]{\n',iDend);
            pos1 = pDend.startPoint*1e3;
            pos2 = pDend.endPoint*1e3;
            fprintf(fid,'pt3dadd(%f,%f,%f,0)\n',pos1(1),pos1(2),pos1(3));
            fprintf(fid,'pt3dadd(%f,%f,%f,0)\n',pos2(1),pos2(2),pos2(3));
            fprintf(fid,'}\n\n');
            
            for iMid = 0:1
                if iMid == 0, mDend = pDend.child1; end
                if iMid == 1, mDend = pDend.child2; end
                fprintf(fid,'middend[%d]{\n',iDend*2+iMid);
                pos1 = mDend.startPoint*1e3;
                pos2 = mDend.endPoint*1e3;
                fprintf(fid,'pt3dadd(%f,%f,%f,0)\n',pos1(1),pos1(2),pos1(3));
                fprintf(fid,'pt3dadd(%f,%f,%f,0)\n',pos2(1),pos2(2),pos2(3));
                fprintf(fid,'}\n\n');
                
                for iDist = 0:1
                    if iDist == 0, dDend = mDend.child1; end
                    if iDist == 1, dDend = mDend.child2; end
                    fprintf(fid,'distdend[%d]{\n',iDend*4+iMid*2+iDist);
                    pos1 = dDend.startPoint*1e3;
                    pos2 = dDend.endPoint*1e3;
                    fprintf(fid,'pt3dadd(%f,%f,%f,0)\n',pos1(1),pos1(2),pos1(3));
                    fprintf(fid,'pt3dadd(%f,%f,%f,0)\n',pos2(1),pos2(2),pos2(3));
                    fprintf(fid,'}\n\n');
                end
            end
        end
    end

    nodectr = 0;
    mysactr = 0;
    flutctr = 0;
    stinctr = 0;
    
    for ctr = 1:lastnode
        ind = axonseg(ctr,7);
        switch ind
            case 1 %node
                fprintf(fid,'node[%d]{\n',nodectr);
                nodectr = nodectr + 1;
            case {2,8} % MYSA
                fprintf(fid,'MYSA[%d]{\n',mysactr);
                mysactr = mysactr + 1;
            case {3,7} % FLUT
                fprintf(fid,'FLUT[%d]{\n',flutctr);
                flutctr = flutctr + 1;
            case {4,5,6}
                fprintf(fid,'STIN[%d]{\n',stinctr);
                stinctr = stinctr + 1;
            otherwise
                error('bad segment index')
        end
        pos1 = axonseg(ctr,:);
        pos2 = axonseg(ctr+1,:);
        fprintf(fid,'pt3dadd(%f,%f,%f,0)\n',pos1(1),pos1(2),pos1(3));
        fprintf(fid,'pt3dadd(%f,%f,%f,0)\n',pos2(1),pos2(2),pos2(3));
        fprintf(fid,'}\n\n');
    end    
    fclose(fid);
    
end
close(h);

end  % function end


% ------------------------------------------------------------------------
% Function to calculate the distance between adjacent points
function [dist]=distcalc(coord1,coord2)
    dist=sqrt((coord2(1)-coord1(1))^2+(coord2(2)-coord1(2))^2+(coord2(3)-coord1(3))^2);
end  % function end


% ------------------------------------------------------------------------
% Function to generate segments from an array of camera lucida points
function [seg] = axonpointparse(tempseg)

    % Define paramter lengths and number of segments/parameter for the axon
    nodeL = 1;      % NODE length
    mysaL = 3;      % MYSA length
    flutL = 10;     % FLUT length
    stinL = 57.67;  % STIN length
    segLength = [nodeL,mysaL,flutL,stinL,stinL,stinL,flutL,mysaL];
    segType = [1,2,3,4,5,6,7,8];    % 0-initial segment 1-node, 2-MYSA, 3-FLUT, 4-STIN...
    
    % Initialize coordinates in 'seg' corresponding to x,y,z
    x1 = tempseg(1,1); x2 = tempseg(2,1);
    y1 = tempseg(1,2); y2 = tempseg(2,2);
    z1 = tempseg(1,3); z2 = tempseg(2,3);
    seg(1,:) = [x1,y1,z1,1,1,0,segType(1)];

    % Define loop variables for parsing through the tree
    k = 1;              % axon points counter
    ctr = 2;            % compartment number counter
    segLctr = 1;        % counter (axon node compartment length)
    segTctr = 2;        % counter (axon node compartment type)
    
    % Parse thru all camera lucida points and create segments of length segL
    while k<size(tempseg,1)-1,

        % Resolve issue if two camera lucida points are at the same coordinates
        while x1==x2 && y1==y2 && z1==z2,
            k = k + 1;
            x1 = tempseg(k,1); x2 = tempseg(k+1,1);
            y1 = tempseg(k,2); y2 = tempseg(k+1,2);
            z1 = tempseg(k,3); z2 = tempseg(k+1,3);
        end

        % Calculate distance between camera lucida points
        tempdist = distcalc([x1 y1 z1],[x2 y2 z2]);
        
        % Increment to next camera lucida point if segment is too short
        while tempdist<segLength(segLctr),
            
            % Resolve the issue of reaching the end of the axon
            if k > length(tempseg)-2, break; end;
            
            % Advance to the next point
            x2 = tempseg(k+1,1);
            y2 = tempseg(k+1,2);
            z2 = tempseg(k+1,3);
            tempdist = distcalc([x1 y1 z1],[x2 y2 z2]);
            k = k + 1;

        end
        
        % Interpolate a segment between coordinates with length segL
        segx = sign(x2-x1)*sqrt(segLength(segLctr)^2*(x2-x1)^2/((x2-x1)^2+(y2-y1)^2+(z2-z1)^2))+x1;
        segy = sign(y2-y1)*sqrt(segLength(segLctr)^2*(y2-y1)^2/((x2-x1)^2+(y2-y1)^2+(z2-z1)^2))+y1;
        segz = sign(z2-z1)*sqrt(segLength(segLctr)^2*(z2-z1)^2/((x2-x1)^2+(y2-y1)^2+(z2-z1)^2))+z1;
        seg(ctr,:) = [segx, segy, segz, 1, ctr, ctr-1, segType(segTctr)];
        x1 = seg(ctr,1); y1 = seg(ctr,2); z1 = seg(ctr,3);

        ctr = ctr + 1;
        
        % If input array is axonal points, determine segment type and define segL accordingly
        segLctr = segLctr + 1; segTctr = segTctr + 1;
        if segLctr > length(segLength),
            segLctr = 1;
        end
        if segTctr > length(segType),
            segTctr = 1;
        end
        
    end

end  % function end
