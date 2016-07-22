function ok=FileExists(name)
    f=fopen(name);
    if f>0
        fclose(f);
        ok=1;
    else
        ok=0;
    end;
