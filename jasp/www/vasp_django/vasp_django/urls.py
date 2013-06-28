from django.conf.urls import patterns, include, url

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

# these have to be in the right order, specific to general
urlpatterns = patterns('',
    url(r'^$',             'vasp.views.index'),
    url(r'^(?P<path>.*)/file/(?P<fname>.*)$', 'vasp.views.get_file'),
    url(r'^(?P<path>.*)/resource/(?P<resource>.*)$', 'vasp.views.get_resource'),
    url(r'^(?P<path>.+)(?P<format>/(html|xml|json|python|text)?/?)$', 'vasp.views.jaspsum'),
    # Examples:
    #url(r'^$', 'vasp_django.views.home', name='home'),
    # url(r'^vasp_django/', include('vasp_django.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    url(r'^admin/', include(admin.site.urls)),
)
