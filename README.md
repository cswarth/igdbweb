# Immunoglobulin Germline Database Website Generator

Generate a static website to navigate the Immunoglobulin Germline
Database (IGDB).  Scripts and templates in this repo build a static
website to help navigate among the Ig germlines.  The raw "database"
is store as JSON formatted objects.  The database is not stored in
this repo - you can find the actual database in the [`IGDB/content`](https://github.com/cswarth/igdb/tree/master/content).

## Usage

To build a static website, assuming your germline database is in
`content/` and your website is hosted out of `/var/www`,

```
    $ pip install git+git://github.com/cswarth/igdbweb@master
    $ python -m igdbweb -c "content/" -o "/var/www"
```

You'll need some raw data to put into `content\`.  The [IGDB repo](https://github.com/cswarth/igdb) holds
the raw germline database.

In production, the process of turning the raw database into a static,
hosted website is run automatically as part of the CI/CD
setup in [IGDB](https://github.com/cswarth/igdb).  See
[`Wercker.yml`](https://github.com/cswarth/igdb/blob/master/wercker.yml)
for the configuration of automatic builds and deployment of the web
site.  
