function ok=DirectoryExists(name)
  q=dir(name);
  ok=(numel(q)>0);
