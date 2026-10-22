#!/usr/bin/env python3

# ======================= python ======================= #
#                         ____ _        ___        __    #
#  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  #
# | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| #
# | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ #
# | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ #
# |_|                                                    #
# ====================================================== #

from dataclasses import dataclass


@dataclass
class CustomSettings:

    FIELD_ORDER = [

    ]

    def values(self):
        return [getattr(self, field) for field in self.FIELD_ORDER]