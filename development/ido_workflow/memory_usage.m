
function [] = memory_usage (mem)
% ANALYZE_MEMORY count memory consumed by workspace variables.
% Input: mem - the command 'whos' from the caller
% GGG how to get an official memory analysis in unix (like "memory" in windows)?

    total_memory = 0;
    for i = 1:size(mem,1)
        total_memory = total_memory + mem(i).bytes;
    end
    
    GB = (2^10)^3;
    log_message('Memory consumed by variables: %.2f GB', total_memory/GB);

end
