function testlog

open_log('loglog.txt');
for k=1:10000
    log_message(sprintf('Message %03d',k));
end
close_log;
