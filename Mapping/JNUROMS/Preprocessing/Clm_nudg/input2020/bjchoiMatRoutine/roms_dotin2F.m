% Script to generate a first draft of the roms code segments that
% read, report, and write to nc files the input parameter list for
% modules such as ecosim

% A simplified version of the ecosim.in file is required for input,
% that has all the preamble and documentation stripped out so it
% consists only of:
% 1-line variable description (that commences with ! in column 1), and 
% 1-line initialization (i.e. the '==' syntax)
% These lines will be parsed to build the code segments
% The output will need editing when inserted inot the roms code, 
% but this does most of the tedious work



% set up the file ids
% stripped down input list
file_in = 'ecosim.in.list';

% outputs
file_inp1 = 'inp_par_segment1';
file_inp2 = 'inp_par_segment2';
file_def = 'def_info_segment';
file_wrt = 'wrt_info_segment';


% open files for writing
fid   = fopen(file_in);
fidi1 = fopen(file_inp1,'w');
fidi2 = fopen(file_inp2,'w');
fidd  = fopen(file_def,'w');
fidw  = fopen(file_wrt,'w');

while 1

  % Parse each variable name and description
  % -----------------------------------------------------------
  
  % read from input list
  % first line should be the description
  
  tline = fgetl(fid);
  % endof file check
  if ~ischar(tline)
    break
  end  
  disp(tline)
  if isempty(tline) 
    % there was a blank line, so read again
    tline = fgetl(fid);
  end
  
  blanks = min(findstr(tline,' '));
  tline(1:blanks) = []; % strip the leading ! and blank spaces
  tline = deblank(tline); % remove trailing blanks
  
  % attempt to parse the units
  bracketsopen = min(findstr(tline,'('));
  bracketsclose = min(findstr(tline,')'));
  if ~isempty(bracketsopen) & ~isempty(bracketsclose)
    varu = tline(bracketsopen+1:bracketsclose-1);
    vard = deblank(tline(1:bracketsopen-1));
  else
    varu = 'unknown';
    vard = tline;
  end
    
  if length(vard)>54
    vard(54:end) = [];
  end
    
  % read from input list
  tline = fgetl(fid);
  if isempty(tline) 
    % there was a blank line, so read again
    tline = fgetl(fid);
  end
  
  % second line should be the variable name (followed by ==
  % initialization)
  
  eq = min(findstr(tline,'='));
  tline(eq:end) = []; % strip everything after the first =
  blanks = findstr(tline,' ');
  tline(blanks) = [];  
  varn = tline;

  disp('...')
  disp([ 'Variable:    '   varn])
  disp([ 'Description: '   vard])
  disp([ 'Units:       '   varu])
  disp(' ')
  
  % code segment for first part in inp_par.F
  % read .in file
  x12  = [ '          '];
  str1 = [x12 'ELSE IF (TRIM(KeyWord).eq.''' varn ''') THEN'];
  str2 = [x12 '  Npts=load_i(Nval, Rval, Ngrids, ' varn ')'];
  
  % code segment for second part of inp_par.F
  % write to stdout
  str3 = [x12 'WRITE (out,100) ' varn '(ng), ''' varn ''','];
  str3 = [str3 repmat(' ',[1 length(length(str3):70)]) '&'];
  x18  = [ '                '];
  str4 = [ x18 '''' vard ''''];
  str4(6) = '&';
  
  % write code segments
  fprintf(fidi1,'%s',str1);
  fprintf(fidi1,'\n');
  fprintf(fidi1,'%s',str2);
  fprintf(fidi1,'\n');
  
  fprintf(fidi2,'%s',str3);
  fprintf(fidi2,'\n');
  fprintf(fidi2,'%s',str4);
  fprintf(fidi2,'\n');

  % code segment for def_info.F
  
  x7 = ['      '];
  str1 = [ x7 'Vinfo( 1)=''' varn ''''];
  str2 = [ x7 'Vinfo( 2)=''' vard ''''];
  str3 = [ x7 'Vinfo( 3)=''' varu ''''];
  str4 = [ x7 'status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)'];
  
  fprintf(fidd,'%s',str1);
  fprintf(fidd,'\n');
  fprintf(fidd,'%s',str2);
  fprintf(fidd,'\n');
  fprintf(fidd,'%s',str3);
  fprintf(fidd,'\n');
  fprintf(fidd,'%s',str4);
  fprintf(fidd,'\n');
  fprintf(fidd,'\n');

  % code segment for wrt_info.F
  
  x9 = ['      '];
  
  str1 = [ x9 'wrt_info=nf_inq_varid(ncid,''' varn '''' ',varid)'];
  str2 = [ x9 'wrt_info=nf_put_var1_TYPE(ncid,varid,1,' varn '(ng))'];
  str3 = [ x9 'IF (wrt_info.ne.nf_noerr) THEN'];
  str4 = [ x9 '  WRITE (stdout,10) ''' varn '''' ', TRIM(ncname)'];
  str5 = [ x9 '  exit_flag=3'];
  str6 = [ x9 '  RETURN'];
  str7 = [ x9 'ENDIF'];

  fprintf(fidw,'%s',str1);
  fprintf(fidw,'\n');
  fprintf(fidw,'%s',str2);
  fprintf(fidw,'\n');
  fprintf(fidw,'%s',str3);
  fprintf(fidw,'\n');
  fprintf(fidw,'%s',str4);
  fprintf(fidw,'\n');
  fprintf(fidw,'%s',str5);
  fprintf(fidw,'\n');
  fprintf(fidw,'%s',str6);
  fprintf(fidw,'\n');
  fprintf(fidw,'%s',str7);
  fprintf(fidw,'\n');
  fprintf(fidw,'\n');

end