function [ ] = DE1(NP,G,D)
    %%%%%%%%%%%%%%差分进化算法求函数极值%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%初始化%%%%%%%%%%%%%%%%%%%%
    MAXRUN=20;
    global fes;
    fes=0;
    ok=0;
    % NP=50;
    % D=10;
    % G=200;
    F0=0.4;
    CR=0.1;
    Xs=20;
    Xx=-20;
    error=10^-6;
    maxFes=100000;
    bestGen=zeros(1,MAXRUN);
    bestFes=zeros(1,MAXRUN);
    bestTime=zeros(1,MAXRUN);
    bestFun=zeros(1,MAXRUN);
    %%%%%%%%%%%%%%%%计算目标函数%%%%%%%%%%%%%%%%%
    if isempty(dir('De1Output'))
    else
        rmdir('De1Output', 's');
    end
    mkdir('De1Output');
    for run = 1:MAXRUN
        tic
        fes=0;
        %%%%%%%%%%%%%%%%%%赋初值%%%%%%%%%%%%%%%%%%%%
        x = zeros(D,NP); %初始种群
        v = zeros(D,NP); %变异种群
        u = zeros(D,NP); %选择种群
        x = rand(D,NP)*(Xs-Xx)+Xx;
        %%%%%%%%%%%%%%%%%开始执行%%%%%%%%%%%%%%%%%%%%
        for m = 1:NP
            Ob(m) = func1(x(:,m));
        end
        trace(1) = min(Ob);
        %%%%%%%%%%%%%%%%%创建文件%%%%%%%%%%%%%%%%%%%%
        fprintf('第%d次运行\t第%d代\t%g\n',run,0,trace(1));
        onerunfile = fopen(['.\De1Output\F1_run' num2str(run) '.txt'],'w'); %建立第run次运行的记录文档
        if fes==NP
            fprintf(onerunfile,'代数\t函数运行次数\t计算最小值\n');
        end
        fprintf(onerunfile,'%-d\t%10d\t\t%g\r\n',0,fes,trace(1));
        %%%%%%%%%%%%%%%差分进化循环%%%%%%%%%%%%%%%%%%
        for gen = 1:G
            %%%%%%%%%%%%%%变异操作%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%r1,r2,r3和m互不相同%%%%%%%%%%%%%%%%%%%%
            lamda = exp(1-G/(G+1-gen));
            F = F0*2^(lamda);
            for m = 1:NP
                r1 = randi([1,NP],1,1);
                while (r1==m)
                    r1 = randi([1,NP],1,1);
                end
                r2 = randi([1,NP],1,1);
                while (r2==m) | (r2==r1)
                    r2 = randi([1,NP],1,1);
                end
                r3 = randi([1,NP],1,1);
                while (r3==m) | (r3==r1) | (r3==r2)
                    r3 = randi([1,NP],1,1);
                end
                v(:,m) = x(:,r1) + F*(x(:,r2)-x(:,r3));
            end

            %%%%%%%%%%%%%%变异算子%%%%%%%%%%%%%%%%%%%%
            r = randi([1,D],1,1);

            for n = 1:D
                cr = rand(1);
                if (cr<=CR) | (n==r)
                    u(n,:) = v(n,:);
                else
                    u(n,:) = x(n,:);
                end        
            end

            %%%%%%%%%%%%%%边界条件的处理%%%%%%%%%%%%%%%%%%%%
            for n = 1:D
                for m = 1:NP
                    if (u(n,m)<Xx) | (u(n,m)>Xs)
                        u(n,m) = rand*(Xs-Xx)+Xx;
                    end
                end
            end

            %%%%%%%%%%%%%%选择操作%%%%%%%%%%%%%%%%%%%%
            for m = 1:NP
                Ob1(m) = func1(u(:,m));
            end
            for m = 1:NP
                if Ob1(m)<Ob(m)
                    x(:,m) = u(:,m);
                end
            end
            for m = 1:NP
                Ob(m) = func1(x(:,m));
            end
            trace(gen+1) = min(Ob);

            if fes >= maxFes
                fprintf(onerunfile,'%d\t%10d\t\t%g\r\n',gen,fes,trace(gen+1));
                break;
            end
            if(trace(gen)~=trace(gen+1) & abs((trace(gen)-trace(gen+1)))<error)
                ok=ok+1;
                break;
            end
            if trace(gen) > trace(gen+1)
                bestGen(run) = gen;
                bestFes(run) = fes;
                bestFun(run) = trace(gen+1);
                bestTime(run) = toc;
            end
            if min(Ob(m)) < error
                break;
            end
            fprintf(onerunfile,'%d\t%10d\t\t%g\r\n',gen,fes,trace(gen+1));
        end
        fclose(onerunfile);
        [SortOb,Index] = sort(Ob);
        x = x(:,Index);
        x = x(:,1);

        Y = min(Ob);
        % figure
        % plot(trace);
        % xlabel('迭代次数');
        % ylabel('目标函数曲线')
        % title('DE目标函数曲线')
    end

    F1total = fopen('.\De1Output\F1total.txt','w');
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
    Total = fopen('.\De1Output\Total.txt','w');
    fprintf(Total,'测试函数名:DE\r\n');
    fprintf(Total,'D:%d\r\n',D);
    fprintf(Total,'Pc:%.1f\r\n',CR);
    fprintf(Total,'Pm:%.1f\r\n',F0);
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

function result = func1(x)
    summ = sum(x.^2);
    result = summ;
    global fes;
    fes = fes+1;
end


