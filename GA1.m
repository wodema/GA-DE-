% function [ ] = GA(NP,G,D)
%     MAXRUN=20;
%     global fes;
%     fes=0;
%     ok=0;
%     error=10^-6;
%     maxFes=100000;
%     bestGen=zeros(1,MAXRUN);
%     bestFes=zeros(1,MAXRUN);
%     bestTime=zeros(1,MAXRUN);
%     bestFun=zeros(1,MAXRUN);    
%     % NP = 50;
%     % D = 10;
%     % G = 1000;
%     Pc = 0.8;
%     Pm = 0.1;
%     Xs = 20;
%     Xx = -20;

%     if isempty(dir('Ga1Output'))
%     else
%         rmdir('Ga1Output', 's');
%     end
%     mkdir('Ga1Output');
%     for run = 1:MAXRUN
%         tic
%         fes=0;
%         %%%%%%%%%%%%%%%%%%����ֵ%%%%%%%%%%%%%%%%%%%%
%         f = zeros(D,NP);
%         nf = zeros(D,NP);
%         f = rand(D,NP)*(Xs-Xx)+Xx;
    
%         for np = 1:NP
%             FIT(np) = func1(f(:,np));
%         end
%         [SortFIT,Index] = sort(FIT);
%         Sortf = f(:,Index);

%         FirstMin = min(FIT);
%         %%%%%%%%%%%%%%%%%�����ļ�%%%%%%%%%%%%%%%%%%%%
%         fprintf('��%d������\t��%d��\t%g\n',run,0,FirstMin);
%         onerunfile = fopen(['.\Ga1Output\F1_run' num2str(run) '.txt'],'w'); %������run�����еļ�¼�ĵ�
%         if fes==NP
%             fprintf(onerunfile,'����\t�������д���\t������Сֵ\n');
%         end
%         fprintf(onerunfile,'%-d\t%10d\t\t%g\r\n',0,fes,FirstMin);
%         for gen = 1:G
%             for i = 1:2:NP
%                 p = rand;
%                 if p < Pc
%                     q = randi([0,1],1,D);
%                     for j = 1:D
%                         if q(j)==1
%                             temp = nf(j,i+1);
%                             nf(j,i+1) = nf(j,i);
%                             nf(j,i) = temp;
%                         end
%                     end 
%                 end
%             end

%             for m = 1:NP
%                 for n = 1:D
%                     r = rand(1,1);
%                     if r < Pm
%                         nf(n,m) = rand(1,1)*(Xs-Xx)+Xx;
%                     end
%                 end
%             end

%             %%%%%%%�������̶ĵĸ��Ʋ���%%%%%%%
%             sum_Fit = sum(SortFIT);
%             fitvalue = SortFIT./sum_Fit;
%             %%%%%%cumsum������ÿ��������ֵ��Ϊ��֮ǰ����ԭ����ֵ֮�͡���%%%%%%
%             fitvalue = cumsum(fitvalue);
%             %%%%%%rand(NP,1)���NP��1����������飬֮������
%             ms = sort(rand(1,NP));
%             fiti = 1;
%             newi = 1;
%             %%%%%%�������临�Ʋ���%%%%%%%
%             while newi <= NP
%                 if(ms(newi)) < fitvalue(fiti)
%                     nf(:,newi) = f(:,fiti);
%                     newi = newi+1;
%                 else
%                     fiti = fiti+1;   
%                 end 
%             end
%             f = nf;
%             f(1,:) = fBest;
%             trace(gen) = min(Fit);


%             if fes >= maxFes
%                 break;
%             end
%             if(gen>1 & trace(gen-1)~=trace(gen) & abs((trace(gen-1)-trace(gen)))<error)
%                 ok=ok+1;
%                 break;
%             end
%             if gen>1 & trace(gen) < trace(gen-1)
%                 bestGen(run) = gen;
%                 bestFes(run) = fes;
%                 bestFun(run) = trace(gen);
%                 bestTime(run) = toc;
%             end
%             fprintf(onerunfile,'%-d\t%10d\t\t%g\r\n',gen,fes,trace(gen));
%         end
%         fclose(onerunfile);
%         % Bestf = Sortf(:,1);
%         % trace(end);
%         % figure;
%         % plot(trace);
%         % xlabel('��������')
%         % ylabel('Ŀ�꺯��ֵ')
%         % title('��Ӧ�Ƚ�������')
%     end
        
%     F1total = fopen('.\Ga1Output\F1total.txt','w');
%     for i = 1:MAXRUN
%         fprintf(F1total,'��%d��:',i);
%         fprintf(F1total,'(�״����Ž�%g)---',bestFun(i));
%         fprintf(F1total,'(������%d)---',bestGen(i));
%         fprintf(F1total,'(�������۴���%d)---',bestFes(i));
%         fprintf(F1total,'(��ʱ%gs)\r\n',bestTime(i));
%     end    
%     fclose(F1total);

%     bestfun=min(bestFun);
%     worstfun=max(bestFun);
%     meanfun=mean(bestFun);
%     stdfun=std(bestFun,0);
%     meanbestGen=mean(bestGen);
%     meanbestFes=mean(bestFes);
%     meanbestTime=mean(bestTime);
%     okPercent=ok/MAXRUN*100;
%     Total = fopen('.\Ga1Output\Total.txt','w');
%     fprintf(Total,'���Ժ�����:DE\r\n');
%     fprintf(Total,'D:%d\r\n',D);
%     fprintf(Total,'Pc:%.1f\r\n',Pc);
%     fprintf(Total,'Pm:%.1f\r\n',Pm);
%     fprintf(Total,'NP:%d\r\n',NP);
%     fprintf(Total,'MAXG:%d\r\n',G);
%     fprintf(Total,'MAXRUN:%d\r\n',MAXRUN);
%     fprintf(Total,'bestfun:%d\r\n',bestfun);
%     fprintf(Total,'worstfun:%d\r\n',worstfun);
%     fprintf(Total,'meanfun:%d\r\n',meanfun);
%     fprintf(Total,'stdfun:%d\r\n',stdfun);
%     fprintf(Total,'meanbestGen:%g\r\n',meanbestGen);
%     fprintf(Total,'meanbestFes:%g\r\n',meanbestFes);
%     fprintf(Total,'meanbestTime:%g\r\n',meanbestTime);
%     fprintf(Total,'�������ٷֱ�(ok):%.2f%%\r\n',okPercent);
%     fclose(Total);
% end


% %%%%%%%%%%��Ӧ�Ⱥ���%%%%%%%%%%
% function result = func1(x)
%     result =  sum(x.^2);
%     global fes;
%     fes = fes+1;
% end

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

    if isempty(dir('Ga1Output'))
    else
        rmdir('Ga1Output', 's');
    end
    mkdir('Ga1Output');
    for run = 1:MAXRUN
        tic
        fes=0;
        f = randi([0,1], NP, D);

        %%%%%%%%%%%%%%%%%�����ļ�%%%%%%%%%%%%%%%%%%%%
        % fprintf('��%d������\t��%d��\t%g\n',run,0,FirstMin);
        onerunfile = fopen(['.\Ga1Output\F1_run' num2str(run) '.txt'],'w'); %������run�����еļ�¼�ĵ�
        % if fes==NP
        %     fprintf(onerunfile,'����\t�������д���\t������Сֵ\n');
        % end
        % fprintf(onerunfile,'%-d\t%10d\t\t%g\r\n',0,fes,FirstMin);
        for gen = 1:G
            for i = 1:NP
                U = f(i,:);
                m = 0;
                for j = 1:D
                    m = U(j)*2^(j-1)+m;
                end
                x(i) = Xx+m*(Xs-Xx)/(2^D-1);
                Fit(i) = func1(x(i));
            end

            maxFit = max(Fit);
            minFit = min(Fit);
            rr = find(Fit == minFit);
            %%%ȡ�׸����Ÿ�������%%%
            fBest = f(rr(1,1),:);
            %%%ȡ�׸�ӳ��������Xֵ%%%
            xBest = x(rr(1,1));
            
            SortFit = sort(Fit);
            %%%%%%%�������̶ĵĸ��Ʋ���%%%%%%%
            sum_Fit = sum(SortFit);
            fitvalue = SortFit./sum_Fit;
            %%%%%%cumsum������ÿ��������ֵ��Ϊ��֮ǰ����ԭ����ֵ֮�͡���%%%%%%
            fitvalue = cumsum(fitvalue);
            %%%%%%rand(NP,1)���NP��1����������飬֮������
            ms = sort(rand(NP,1));
            fiti = 1;
            newi = 1;
            %%%%%%�������临�Ʋ���%%%%%%%
            while newi <= NP
                if(ms(newi)) < fitvalue(fiti)
                    nf(newi,:) = f(fiti,:);
                    newi = newi+1;
                else
                    fiti = fiti+1;   
                end 
            end
        
            %%%%%%%�������%%%%%%%%
            for i = 1:2:NP
                p = rand;
                if p < Pc
                    q = randi([0,1],1,D);
                    for j = 1:D
                        if q(j)==1
                            temp = nf(i+1,j);
                            nf(i+1,j) = nf(i,j);
                            nf(i,j) = temp;
                        end
                    end 
                end
            end
        
            %%%%%%%�������%%%%%%%
            i = 1;
            while i <= round(NP*Pm)
                h = randi(NP);
                for j = 1:round(D*Pm)
                    g = randi(D);
                    nf(h,g) =~ nf(h,g);
                end
                i = i+1;
            end
            f = nf;
            f(1,:) = fBest;
            trace(gen) = minFit;

            if fes >= maxFes
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
        % xlabel('��������')
        % ylabel('Ŀ�꺯��ֵ')
        % title('��Ӧ�Ƚ�������')
    end
    figure;
    plot(trace); 
    xlabel('��������')
    ylabel('Ŀ�꺯��ֵ')
    title('��Ӧ�Ƚ�������')

    F1total = fopen('.\Ga1Output\F1total.txt','w');
    for i = 1:MAXRUN
        fprintf(F1total,'��%d��:',i);
        fprintf(F1total,'(�״����Ž�%g)---',bestFun(i));
        fprintf(F1total,'(������%d)---',bestGen(i));
        fprintf(F1total,'(�������۴���%d)---',bestFes(i));
        fprintf(F1total,'(��ʱ%gs)\r\n',bestTime(i));
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
    Total = fopen('.\Ga1Output\Total.txt','w');
    fprintf(Total,'���Ժ�����:DE\r\n');
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
    fprintf(Total,'�������ٷֱ�(ok):%.2f%%\r\n',okPercent);
    fclose(Total);
end


%%%%%%%%%%��Ӧ�Ⱥ���%%%%%%%%%%
function result = func1(x)
    result =  sum(x.^2);
    global fes;
    fes = fes+1;
end