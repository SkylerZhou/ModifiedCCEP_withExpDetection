% implements an additional constriant where both N1 and N2 components of
% the CCEPs must be present and valid for a signal to be kept. 
function out = AA_require_both_Ns(out)

NumRows = size(out.rejection_details(1).reject.keep,1);
NumColumns = size(out.rejection_details(1).reject.keep,2);

for ich=1:NumRows
    ich;
    for jch=1:NumColumns
        if isempty(out.elecs(ich).arts)
            out.rejection_details(1).reject.keep(ich,jch) = 0;
            out.rejection_details(2).reject.keep(ich,jch) = 0;
            out.rejection_details(1).reject.empty(ich,jch)= 1;
            out.rejection_details(2).reject.empty(ich,jch) = 1;
        end
        if size(out.elecs(ich).N1,1) >=2
            if isnan(out.elecs(ich).N1(jch,1)) || isnan(out.elecs(ich).N2(jch,1)) || (out.elecs(ich).N1(jch,1)) == 0 || (out.elecs(ich).N2(jch,1)) == 0
                out.rejection_details(1).reject.keep(ich,jch) = 0;
                out.rejection_details(2).reject.keep(ich,jch) = 0;
                out.rejection_details(1).reject.no_both(ich,jch)= 1;
                out.rejection_details(2).reject.no_both(ich,jch) = 1;
            end
            % if size(out.elecs(ich).detrend_filt_avgs,1)<=2
            %     out.rejection_details(1).reject.empty(ich,jch)=1;
            %     out.rejection_details(2).reject.empty(ich,jch)=1;
            %     out.rejection_details(1).reject.keep(ich,jch) = 0;
            %     out.rejection_details(2).reject.keep(ich,jch) = 0;
            % end
            % if out.rejection_details(1).reject.sig_avg(ich,jch) == 1 || out.rejection_details(1).reject.pre_thresh(ich,jch) == 1 || out.rejection_details(1).reject.at_thresh(ich,jch) ==1
            %     out.rejection_details(1).reject.no_both(ich,jch) == 0;
            % end
        end
        % if size(out.elecs(ich).N1,1) <=2
        %     out.rejection_details(1).reject.keep(ich,jch) = 0;
        %     out.rejection_details(2).reject.keep(ich,jch) = 0;
        %     out.rejection_details(1).reject.empty(ich,jch)=1;
        %     out.rejection_details(2).reject.empty(ich,jch)=1;
        % end
        % if out.rejection_details(1).reject.keep(ich,jch) ==0 || out.rejection_details(2).reject.keep(ich,jch) ==0
        %     out.rejection_details(1).reject.keep(ich,jch) = 0;
        %     out.rejection_details(2).reject.keep(ich,jch) = 0;
        %     out.rejection_details(1).reject.empty(ich,jch)=1;
        %     out.rejection_details(2).reject.empty(ich,jch)=1;
        % end
    end
end
end