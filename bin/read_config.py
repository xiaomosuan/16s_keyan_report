import configparser


class MyConfigParser(configparser.ConfigParser):
    """
    自定义MyConfigParser， 继承configparser.ConfigParser， 解决option被自动转为小写
    """
    def init(self, defaults=None):
        configparser.ConfigParser.init(self, defaults=defaults)

    def optionxform(self, optionstr):
        return optionstr


def get_config(config_file):
    config = MyConfigParser()
    config.read(config_file)
    return config


def get_items(sec, config_file='/Users/congliu/OneDrive/普瑞森/New_keyan_report/ini/config_XIAOJIBASHI.ini'):
    config = get_config(config_file)
    items = config.items(sec)
    return items

