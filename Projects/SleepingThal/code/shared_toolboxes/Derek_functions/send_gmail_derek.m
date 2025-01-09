function [] = send_gmail_derek(subject, message)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    fpath_to_gmail_passcode = "Y:\matlab_gmail_passcode.txt";
    password = string(fileread(fpath_to_gmail_passcode));

    mail = 'dossneurosurgicalassociates@gmail.com'; % my gmail address
    host = 'smtp.gmail.com';
    sendto = 'derekjdoss@gmail.com';
    % sendto = 'derek.j.doss@vanderbilt.edu';
    Subject = subject;
    Message = message;

    % preferences
    setpref('Internet','SMTP_Server', host);
    setpref('Internet','E_mail', mail);
    setpref('Internet','SMTP_Username',mail);
    setpref('Internet','SMTP_Password',password);
    props = java.lang.System.getProperties;
    props.setProperty('mail.smtp.auth','true');
    props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
    props.setProperty('mail.smtp.socketFactory.port','465');
    
    % execute
    sendmail(sendto,Subject,Message)

end