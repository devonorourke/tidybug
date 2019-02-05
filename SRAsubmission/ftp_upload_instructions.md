installed FTP with brew.
went to SRA upload site and selected `FTP upload` which had pre-generated tmp username/password info
then used Terminal to connect to FTP site, starting at directory where files were downloaded

```
ftp -i ftp-private.ncbi.nlm.nih.gov   ## logs into FTP site
subftp    ## username
w4pYB9VQ  ## password
cd uploads/devon.orourke_gmail.com_SFmljtHQ   ## go to subfolder
mkdir libA    ## make subfolder within this subfolder
cd libA     ## move into subdirectory
passive on  ## avoids firewall issues
mput *.gz    ## upload all the .gz files (don't just use 'put')
```
This will work; just takes forever.

---

Alternatively, use aspera:
```
/Users/do/Applications/Aspera_connect/Contents/Resources/ascp \
-i /Users/do/.ssh/id_rsa -QT -l100m -k1 -d \
/Users/do/Desktop/tempdata \
subasp@upload.ncbi.nlm.nih.gov:uploads/devon.orourke_gmail.com_k7RfwjBu
```
this command is failing to authenticate though...
