function [ ] = GA(NP,G,D)
    MAXRUN=20;
    global fes;
    fes=0;
    ok=0;
    error=10^-6;
    maxFes=100000;
    bestGen=zeros(1,MAXRUN);
    bestFes=zeros(1,MAXRUN);
    bestTime=zeros(1,MAXRUN);
    bestFun=zeros(1,MAXRUN);    
    % NP = 50;
    % D = 10;
    % G = 1000;
    Pc = 0.8;
    Pm = 0.1;
    Xs = 20;
    Xx = -20;

    if isempty(dir('Ga2Output'))
    else
        rmdir('Ga2Output', 's');
    end
    mkdir('Ga2Output');
    for run = 1:MAXRUN
        tic
        fes=0;
        %%%%%%%%%%%%%%%%%%赋初值%%%%%%%%%%%%%%%%%%%%
        f = zeros(D,NP);
        nf = zeros(D,NP);
        f = rand(D,NP)*(Xs-Xx)+Xx;
    
        for np = 1:NP
            FIT(np) = func1(f(:,np));
        end
        [SortFIT,Index] = sort(FIT);
        Sortf = f(:,Index);

        FirstMin = min(FIT);
        %%%%%%%%%%%%%%%%%创建文件%%%%%%%%%%%%%%%%%%%%
        fprintf('第%d次运行\t第%d代\t%g\n',run,0,FirstMin);
        onerunfile = fopen(['.\Ga2Output\F1_run' num2str(run) '.txt'],'w'); %建立第run次运行的记录文档
        if fes==NP
            fprintf(onerunfile,'代数\t函数运行次数\t计算最小值\n');
        end
        fprintf(onerunfile,'%-d\t%10d\t\t%g\r\n',0,fes,FirstMin);
        for gen = 1:G
            Emper = Sortf(:,1);
            NoPoint = round(D*Pc);
            PoPoint = randi([1 D],NoPoint,NP/2);
            nf = Sortf;
            for i = 1:NP/2
                nf(:,2*i-1) = Emper;
                nf(:,2*i) = Sortf(:,2*i);
                for k = 1:NoPoint
                    nf(PoPoint(k,i), 2*i-1) = nf(PoPoint(k,i),2*i);
                    nf(PoPoint(k,i), 2*i) = Emper(PoPoint(k,i));
                end
            end
    
            for m = 1:NP
                for n = 1:D
                    r = rand(1,1);
                    if r < Pm
                        nf(n,m) = rand(1,1)*(Xs-Xx)+Xx;
                    end
                end
            end
    
            for np = 1:NP
                NFIT(np) = func1(nf(:,np));
            end
            [NSortFIT,Index] = sort(NFIT);
            NSortf = nf(:,Index);
    
            f1 = [Sortf,NSortf];
            FIT1 = [SortFIT,NSortFIT];
            [SortFIT1, Index] = sort(FIT1);
            Sortf1 = f1(:,Index);
            SortFIT = SortFIT1(1:NP);
            Sortf = Sortf1(:,1:NP);
            trace(gen) = SortFIT(1);

            if fes >= maxFes
                fprintf(onerunfile,'%-d\t%10d\t\t%g\r\n',gen,fes,trace(gen));
                break;
            end
            if(gen>1 & trace(gen-1)~=trace(gen) & abs((trace(gen-1)-trace(gen)))<error)
                ok=ok+1;
                break;
            end
            if gen>1 & trace(gen) < trace(gen-1)
                bestGen(run) = gen;
                bestFes(run) = fes;
                bestFun(run) = trace(gen);
                bestTime(run) = toc;
            end

            fprintf(onerunfile,'%-d\t%10d\t\t%g\r\n',gen,fes,trace(gen));
        end
        fclose(onerunfile);
        % Bestf = Sortf(:,1);
        % trace(end);
        % figure;
        % plot(trace);
        % xlabel('迭代次数')
        % ylabel('目标函数值')
        % title('适应度进化曲线')
    end
        
    F1total = fopen('.\Ga2Output\F1total.txt','w');
    for i = 1:MAXRUN
        fprintf(F1total,'第%d代:',i);
        fprintf(F1total,'(首次最优解%g)---',bestFun(i));
        fprintf(F1total,'(迭代数%d)---',bestGen(i));
        fprintf(F1total,'(函数评价次数%d)---',bestFes(i));
        fprintf(F1total,'(用时%gs)\r\n',bestTime(i));
    end    
    fclose(F1total);

    bestfun=min(bestFun);
    worstfun=max(bestFun);
    meanfun=mean(bestFun);
    stdfun=std(bestFun,0);
    meanbestGen=mean(bestGen);
    meanbestFes=mean(bestFes);
    meanbestTime=mean(bestTime);
    okPercent=ok/MAXRUN*100;
    Total = fopen('.\Ga2Output\Total.txt','w');
    fprintf(Total,'测试函数名:DE\r\n');
    fprintf(Total,'D:%d\r\n',D);
    fprintf(Total,'Pc:%.1f\r\n',Pc);
    fprintf(Total,'Pm:%.1f\r\n',Pm);
    fprintf(Total,'NP:%d\r\n',NP);
    fprintf(Total,'MAXG:%d\r\n',G);
    fprintf(Total,'MAXRUN:%d\r\n',MAXRUN);
    fprintf(Total,'bestfun:%d\r\n',bestfun);
    fprintf(Total,'worstfun:%d\r\n',worstfun);
    fprintf(Total,'meanfun:%d\r\n',meanfun);
    fprintf(Total,'stdfun:%d\r\n',stdfun);
    fprintf(Total,'meanbestGen:%g\r\n',meanbestGen);
    fprintf(Total,'meanbestFes:%g\r\n',meanbestFes);
    fprintf(Total,'meanbestTime:%g\r\n',meanbestTime);
    fprintf(Total,'符合误差百分比(ok):%.2f%%\r\n',okPercent);
    fclose(Total);
end


%%%%%%%%%%适应度函数%%%%%%%%%%
function result = func1(x)
    result =  sum(x.^2);
    global fes;
    fes = fes+1;
end