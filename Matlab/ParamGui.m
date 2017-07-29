function ParamGui
format compact
close all
myGui = gui.autogui('Location','float');
myGui.PanelWidth=200;

figure; MeshAxes=axes; axis equal
figure; ParamAxes=axes; axis equal

M=TriangleMesh;
K=TriangleMesh;

axes(MeshAxes);
if verLessThan('matlab','8.4')
    Mh=M.draw('HitTest','off','CDataMapping','scaled');
else
    Mh=M.draw('HitTest','off','CDataMapping','scaled','PickableParts','None');
end
axes(ParamAxes);
Kh=K.draw('HitTest','off','CDataMapping','scaled','PickableParts','None');

error=[];
Parameterizer = ParameterizerClass;
FileName = [];
%% GUI initialization
TxtMenuEnergyExpType = gui.textmenu('Energy type',Parameterizer.EnergyTypesList);
TxtMenuEnergyExpType.Value=Parameterizer.EnergyTypesList{1};
TxtMenuEnergyExpType.ValueChangedFcn = @OnUpdateExpType;

BtnLoadSource = gui.pushbutton('Load Mesh');
BtnLoadSource.ValueChangedFcn = @OnButtonLoadMesh;

BtnInitialize = gui.pushbutton('Initialize');
BtnInitialize.ValueChangedFcn = @OnBtnInitialize;

BtnInitializeFromFile = gui.pushbutton('Initialize From File');
BtnInitializeFromFile.ValueChangedFcn = @OnBtnInitializeFromFile;

TxtMenuProxyType = gui.textmenu('Hessian type',Parameterizer.ProxyTypesList);
TxtMenuProxyType.Value=Parameterizer.ProxyTypesList{1};
TxtMenuProxyType.ValueChangedFcn = @OnUpdateParams;

BtnSolveIter = gui.pushbutton('Do Iteration');
BtnSolveIter.ValueChangedFcn = @OnBtnDoIteration;

BtnSolve = gui.pushbutton('RunExperiment');
BtnSolve.ValueChangedFcn = @OnRunExperiment;

BtnSolveAll = gui.pushbutton('RunExperimentAll');
BtnSolveAll.ValueChangedFcn = @OnRunExperimentAll;

CheckDrawExp = gui.checkbox('Draw Exp');
CheckDrawExp.Value = true;

%% Callbacks
    function OnUpdateExpType(~)
        Parameterizer.SetEnergyExp(TxtMenuEnergyExpType.Value);
        %%TxtMenuProxyType -change this
        TxtMenuProxyType.MenuItems = Parameterizer.ProxyTypesList;
        if strcmp(TxtMenuEnergyExpType.Value,'Neo-Hookean')
           OnButtonLoadMesh;
        end
    end

    function OnButtonLoadMesh(~)
        switch TxtMenuEnergyExpType.Value
            case 'Sym-Dirichlet'
                [MeshFile,PathName] = uigetfile({'*.obj;*.off;*.png'});
                FileName = [PathName,MeshFile];
                [~,~,ext] = fileparts(MeshFile);
                M=TriangleMesh(ext(2:end),FileName);               
            case 'Neo-Hookean'
                M = TriangleMesh('off','bar.off');
            otherwise
                error('incorrect energy type')
        end
            M.centralize('BoundingBox');
            M.V=M.V/sqrt(M.computeSurfaceArea);        
        K = Parameterizer.SetMesh(M); % already initializes
        Draw;
    end

    function OnBtnInitialize(~)
        K = Parameterizer.Initialize;
        Draw;
    end

    function OnBtnInitializeFromFile(~)
        [MeshFile,PathName] = uigetfile({'*.obj;*.off;*.mat'});
        [~,~,ext] = fileparts(MeshFile);
        switch ext
            case '.obj'
                temp=TriangleMesh(ext(2:end),[PathName,MeshFile],'squeeze');
                K.V=temp.Material.Vt;
            case '.off'
                temp=TriangleMesh(ext(2:end),[PathName,MeshFile],'squeeze');
                K.V=temp.V(1:2,:);
            case '.mat'
                S = load([PathName,MeshFile]);
                K.V=reshape(S.x,[],2)';
        end
        Parameterizer.ResetOptimization;
        Draw;
    end

    function OnUpdateParams(~)
        Parameterizer.SetProxyType(TxtMenuProxyType.Value);
    end
    function OnBtnDoIteration(~)
        K.V = Parameterizer.DoIteration(K.V')';
        Draw
    end

    function OnRunExperiment(~)
        RunExperiment({TxtMenuProxyType.Value});
    end

    function OnRunExperimentAll(~)
        ProxyTypesList=Parameterizer.ProxyTypesList;
        RunExperiment(ProxyTypesList);
    end
%% Helper functions
    function RunExperiment(ProxyTypesList)
        n=150;
        data=struct('type',{},'Ki',{},'f',{},'stepSize',{},'g',{},'dx',{},'status',{},'conv_iter',{},'tsession',{});
        [K0,f0,g0] = Parameterizer.Initialize;
        
        for i=1:length(ProxyTypesList)
            tic
            display('________________________________________________________________________________________________________________');
            proxytype=ProxyTypesList{i}
            Ki=zeros(numel(K.V),n);
            f=zeros(1,n);
            stepSize=zeros(1,n);
            f(1)=f0;
            stepSize(1)=0;
            g=zeros(numel(K.V),n);
            status = cell(1,n);
            g(:,1)=g0;
            dx=zeros(1,n);
            Ki(:,1)=reshape(K0.V',[],1);
            Parameterizer.SetProxyType(proxytype);
            Parameterizer.ResetOptimization;
            conv_iter = n;
            bConverged = false;
            for j=2:n
                try
                    [Ki(:,j),f(j),g(:,j),dx(j),stepSize(j)]=Parameterizer.DoIteration(Ki(:,j-1));
                    if CheckDrawExp.Value
                        K.draw(Kh);
                        drawnow;
                    end
                catch
                    conv_iter = j;
                    status{j} = 'Error';
                    break;
                end
                fcur = f(j);
                fprev = f(j-1);
                xcur = Ki(:,j);
                xprev = Ki(:,j-1);
                
                status{j} = Parameterizer.OptimizationConverged(fcur , fprev, xcur, xprev);
                if (~bConverged) && (fprev - fprev < 0)
                    conv_iter = j;
                    status{j} = 'ObjIncreased';
                    bConverged = true;
                    break;
                end
                if (~bConverged) && (~strcmp(status{j},'NotConverged'))
                    conv_iter = j;
                    bConverged = true;
                    break;
                end
            end
            t = toc;
            data{i}.type = proxytype;
            data{i}.Ki = Ki;
            data{i}.f = f(1:conv_iter);
            data{i}.stepSize = stepSize;
            data{i}.g = g;
            data{i}.dx = dx;
            data{i}.status = status;
            data{i}.conv_iter = conv_iter;
            data{i}.tsession = t;
            
            save('curr_data.mat','data');
            
            display(['session time: ' num2str(t)])
            display(['convergencr iter: ' num2str(conv_iter)])
            display(['objective val: ' num2str(f(conv_iter))])
            display(status{conv_iter})
            
        end
        drawGraph(data);
        
    end


    function Draw
        M.draw(Mh);
        K.draw(Kh);
        drawnow;
    end

    function drawGraph(data)
        methods_num = length(data);
        
        %% experiment structure
        % data
        
        %% methods comparison
        expNames = {};
        Nv = num2str(length(data{1}.Ki)/2);
        fig_E_Iter = figure('Name',['Energy - Iter, Nv = ',Nv]);
        ax = gca;
        
        minf = inf;
        for ii=1:methods_num
            if minf>data{ii}.f(data{ii}.conv_iter)
                minf =  data{ii}.f(data{ii}.conv_iter);
            end
        end
        
        for ii=1:methods_num
            exp = data{ii};
            
            % convergence
            display([exp.type,', convergence: ',exp.status{exp.conv_iter} ', num iter: ' num2str(exp.conv_iter) ',  object val: ' num2str(exp.f(exp.conv_iter))]);
            
            f = exp.f-minf+1e-6;
            
            figure(fig_E_Iter)
            hold on;
            
            semilogy(ax,0:exp.conv_iter-1,f(1:exp.conv_iter), 'LineWidth', 1);
            set(gca,'yscale','log')
            
            %             plot(1:exp.conv_iter,exp.f(1:exp.conv_iter),'-');
            
            
            expNames{end+1} = exp.type;
            
        end
        
        figure(fig_E_Iter)
        xlabel('iterations')
        ylabel('Etot')
        legend(expNames{:})
        title('Etot(it)')
        
    end
end